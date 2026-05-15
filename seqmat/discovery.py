"""Gene discovery: list, count, search, summarize the genes installed for an organism.

All queries hit ``genes.db`` directly (fast, accurate), so they work on modern
SQLite-backed installs. Functions return empty containers when an organism is not
configured rather than raising — they're meant for use in notebooks and CLIs where
soft failure reads better than a stack trace.
"""
from __future__ import annotations

import logging
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional

from .config import get_available_organisms, get_organism_config, load_config
from .sqlite_store import get_genes_db_path

_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Organism listings
# ---------------------------------------------------------------------------

def list_available_organisms() -> List[str]:
    """Organism IDs that are currently configured/installed."""
    return get_available_organisms()


def list_supported_organisms() -> List[str]:
    """Alias of :func:`list_available_organisms` (kept for back-compat)."""
    return get_available_organisms()


# ---------------------------------------------------------------------------
# genes.db query helpers
# ---------------------------------------------------------------------------

def _connect(organism: Optional[str]) -> Optional[sqlite3.Connection]:
    """Open a read-only connection to genes.db, or return None if unavailable."""
    db_path = get_genes_db_path(organism)
    if not db_path:
        return None
    try:
        return sqlite3.connect(db_path)
    except sqlite3.Error as exc:
        _log.warning("Failed to open %s: %s", db_path, exc)
        return None


def list_gene_biotypes(organism: Optional[str] = None) -> List[str]:
    """List distinct ``gene_biotype`` values present in ``genes.db``."""
    conn = _connect(organism)
    if conn is None:
        return []
    try:
        rows = conn.execute("SELECT DISTINCT biotype FROM genes ORDER BY biotype").fetchall()
    finally:
        conn.close()
    return [r[0] for r in rows if r[0]]


def count_genes(organism: Optional[str] = None, biotype: Optional[str] = None) -> Dict[str, int]:
    """Count genes per biotype.

    If ``biotype`` is given, returns ``{biotype: count}`` for that one biotype.
    Otherwise returns a dict of all biotypes -> counts.
    """
    conn = _connect(organism)
    if conn is None:
        return {}
    try:
        if biotype:
            row = conn.execute(
                "SELECT COUNT(*) FROM genes WHERE biotype = ?", (biotype,)
            ).fetchone()
            return {biotype: int(row[0])} if row else {biotype: 0}
        rows = conn.execute(
            "SELECT biotype, COUNT(*) FROM genes GROUP BY biotype ORDER BY biotype"
        ).fetchall()
    finally:
        conn.close()
    return {bt: int(n) for bt, n in rows if bt}


def get_gene_list(
    organism: Optional[str] = None,
    biotype: Optional[str] = None,
    limit: Optional[int] = None,
) -> List[str]:
    """List gene names (skipping empty Ensembl symbols), optionally filtered by biotype."""
    conn = _connect(organism)
    if conn is None:
        return []
    sql = "SELECT gene_name FROM genes WHERE gene_name != ''"
    params: tuple = ()
    if biotype:
        sql += " AND biotype = ?"
        params = (biotype,)
    sql += " ORDER BY gene_name"
    if limit:
        sql += " LIMIT ?"
        params = (*params, int(limit))
    try:
        rows = conn.execute(sql, params).fetchall()
    finally:
        conn.close()
    return [r[0] for r in rows]


def get_all_genes(
    organism: Optional[str] = None,
    biotype: Optional[str] = None,
) -> List[Dict[str, str]]:
    """Return all genes for an organism as ``[{organism, biotype, gene_name, gene_id}, ...]``."""
    organism = organism or "<default>"
    conn = _connect(organism if organism != "<default>" else None)
    if conn is None:
        return []
    sql = "SELECT gene_name, gene_id, biotype FROM genes"
    params: tuple = ()
    if biotype:
        sql += " WHERE biotype = ?"
        params = (biotype,)
    sql += " ORDER BY gene_name"
    try:
        rows = conn.execute(sql, params).fetchall()
    finally:
        conn.close()
    return [
        {"organism": organism, "biotype": bt, "gene_name": gn, "gene_id": gid}
        for gn, gid, bt in rows
    ]


def search_genes(
    organism: Optional[str] = None,
    query: str = "",
    biotype: Optional[str] = None,
    limit: int = 10,
) -> List[Dict[str, str]]:
    """Search genes by ``gene_name`` or ``gene_id`` substring (case-insensitive)."""
    if not query:
        return []
    conn = _connect(organism)
    if conn is None:
        return []
    pattern = f"%{query}%"
    sql = (
        "SELECT gene_name, gene_id, biotype FROM genes "
        "WHERE (UPPER(gene_name) LIKE UPPER(?) OR UPPER(gene_id) LIKE UPPER(?))"
    )
    params: tuple = (pattern, pattern)
    if biotype:
        sql += " AND biotype = ?"
        params = (*params, biotype)
    sql += " ORDER BY gene_name LIMIT ?"
    params = (*params, int(limit))
    try:
        rows = conn.execute(sql, params).fetchall()
    finally:
        conn.close()
    org_label = organism or "<default>"
    return [
        {"organism": org_label, "biotype": bt, "gene_name": gn, "gene_id": gid}
        for gn, gid, bt in rows
    ]


def available_genes(organism: str = "hg38", biotype: Optional[str] = "protein_coding") -> Iterator[str]:
    """Yield distinct gene symbols installed for an organism.

    Defaults to protein-coding genes (matches the legacy filesystem-backed behavior).
    Pass ``biotype=None`` to stream every biotype. Empty Ensembl symbols are skipped.
    """
    conn = _connect(organism)
    if conn is None:
        return
    sql = "SELECT DISTINCT gene_name FROM genes WHERE gene_name != ''"
    params: tuple = ()
    if biotype is not None:
        sql += " AND biotype = ?"
        params = (biotype,)
    sql += " ORDER BY gene_name"
    try:
        for (name,) in conn.execute(sql, params):
            yield name
    finally:
        conn.close()


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def get_organism_info(organism: str) -> Dict[str, Any]:
    """Detailed information about an organism: paths + biotype-keyed gene counts."""
    try:
        config = get_organism_config(organism)
    except ValueError:
        return {"error": f"Organism '{organism}' not configured"}

    info: Dict[str, Any] = {
        "organism": organism,
        "configured": True,
        "paths": {k: str(v) for k, v in config.items()},
        "data_available": {},
    }

    counts = count_genes(organism)
    if counts:
        info["data_available"]["biotypes"] = sorted(counts.keys())
        info["data_available"]["gene_counts"] = counts

    chrom_path = config.get("CHROM_SOURCE")
    if chrom_path and Path(chrom_path).exists():
        chrom_files = list(Path(chrom_path).glob("*.fasta"))
        if chrom_files:
            info["data_available"]["chromosomes"] = [f.stem for f in chrom_files]

    return info


def data_summary() -> Dict[str, Any]:
    """A complete data overview across every configured organism."""
    configured = list_available_organisms()
    summary: Dict[str, Any] = {
        "supported_organisms": configured,
        "configured_organisms": configured,
        "organisms": {},
    }
    for organism in configured:
        try:
            summary["organisms"][organism] = get_organism_info(organism)
        except Exception as e:  # pragma: no cover - defensive
            summary["organisms"][organism] = {"error": f"Configuration error: {e}"}

    total_genes = 0
    total_biotypes: set = set()
    for org_info in summary["organisms"].values():
        gc = org_info.get("data_available", {}).get("gene_counts", {})
        for biotype, count in gc.items():
            total_genes += count
            total_biotypes.add(biotype)
    summary["totals"] = {
        "organisms": len(configured),
        "biotypes": len(total_biotypes),
        "genes": total_genes,
    }
    return summary


def print_data_summary() -> None:
    """Print a formatted summary of all installed genomics data."""
    from .config import DEFAULT_ORGANISM_DATA

    summary = data_summary()
    print("SeqMat Genomics Data Summary")
    print("=" * 40)
    totals = summary["totals"]
    print(f"Total: {totals['organisms']} organism(s), {totals['biotypes']} biotype(s), {totals['genes']:,} gene(s)\n")

    print("Configured organisms:")
    configured = set(summary["configured_organisms"])
    for org in sorted(DEFAULT_ORGANISM_DATA.keys() | configured):
        status = "configured" if org in configured else "not configured"
        name = DEFAULT_ORGANISM_DATA.get(org, {}).get("name", org)
        print(f"  - {org}: {name}  [{status}]")
    print()

    for organism, info in summary["organisms"].items():
        print(f"{organism.upper()} data:")
        if "error" in info:
            print(f"  error: {info['error']}")
            continue
        data_avail = info.get("data_available", {})
        if "gene_counts" in data_avail:
            print("  Gene types:")
            for biotype, count in sorted(data_avail["gene_counts"].items()):
                print(f"    {biotype}: {count:,}")
        if "chromosomes" in data_avail:
            chroms = data_avail["chromosomes"]
            preview = ", ".join(sorted(chroms)[:5]) + ("..." if len(chroms) > 5 else "")
            print(f"  Chromosomes ({len(chroms)}): {preview}")
        print("  Paths:")
        for path_name, path_value in info["paths"].items():
            ok = "✓" if Path(path_value).exists() else "✗"
            print(f"    {ok} {path_name}: {path_value}")
        print()


__all__ = [
    "list_available_organisms",
    "list_supported_organisms",
    "list_gene_biotypes",
    "count_genes",
    "get_gene_list",
    "get_all_genes",
    "search_genes",
    "available_genes",
    "get_organism_info",
    "data_summary",
    "print_data_summary",
]
