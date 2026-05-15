"""Fast (chromosome, position) -> gene lookup.

Strategy: build a compact per-chromosome interval index from genes.db, persist it as a
sidecar ``gene_locations.npz`` next to genes.db, and cache it in memory. Queries are
O(log n) via ``np.searchsorted`` plus a vectorized end-bound check on the surviving
candidates — sub-microsecond per lookup for typical chromosomes.

The sidecar is built lazily on first query if missing, and is also produced by
``retrieve_and_parse_ensembl_annotations`` for fresh builds.
"""
from __future__ import annotations

import pickle
import sqlite3
import threading
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

from .config import get_default_organism, get_organism_config


PosArg = Union[int, Tuple[int, int], List[int]]


_INDEX_CACHE: Dict[str, "_LocationIndex"] = {}
_BUILD_LOCK = threading.Lock()


def _normalize_chrm(chrm: str) -> str:
    """Strip a leading 'chr'/'CHR' and uppercase. genes.db stores chromosomes without 'chr'."""
    c = str(chrm).strip()
    if c[:3].lower() == "chr":
        c = c[3:]
    return c.upper()


def _coerce_pos(pos: PosArg) -> Tuple[int, int]:
    """Return (start, end) inclusive bounds for a point or range query."""
    if isinstance(pos, (tuple, list)):
        if len(pos) != 2:
            raise ValueError(f"Range query must be (start, end); got {pos!r}")
        a, b = int(pos[0]), int(pos[1])
        return (a, b) if a <= b else (b, a)
    p = int(pos)
    return p, p


def _locations_path(organism: str) -> Path:
    """Sidecar path: ``<organism dir>/gene_locations.npz``."""
    config = get_organism_config(organism)
    db_path = Path(config["genes_db"])
    return db_path.parent / "gene_locations.npz"


def _genes_db_path(organism: str) -> Path:
    config = get_organism_config(organism)
    return Path(config["genes_db"])


class _LocationIndex:
    """Per-organism index: chrm -> (starts, ends, names) sorted by start."""

    __slots__ = ("by_chrm",)

    def __init__(self, by_chrm: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]]):
        self.by_chrm = by_chrm

    def query_names(self, chrm: str, pos: PosArg) -> List[str]:
        """Return gene names whose interval overlaps the point or range."""
        key = _normalize_chrm(chrm)
        entry = self.by_chrm.get(key)
        if entry is None:
            return []
        starts, ends, names = entry
        qs, qe = _coerce_pos(pos)
        # Candidates: starts <= qe. searchsorted returns the insertion point on the right
        # so all indices [0, hi) satisfy starts <= qe.
        hi = int(np.searchsorted(starts, qe, side="right"))
        if hi == 0:
            return []
        mask = ends[:hi] >= qs
        if not mask.any():
            return []
        return names[:hi][mask].tolist()

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        chroms = sorted(self.by_chrm.keys())
        out: Dict[str, np.ndarray] = {"__chroms__": np.array(chroms, dtype=object)}
        for c in chroms:
            starts, ends, names = self.by_chrm[c]
            out[f"s::{c}"] = starts
            out[f"e::{c}"] = ends
            out[f"n::{c}"] = names
        # allow_pickle is needed because gene-name arrays are object dtype.
        np.savez(str(path), **out)

    @classmethod
    def load(cls, path: Path) -> "_LocationIndex":
        with np.load(str(path), allow_pickle=True) as z:
            chroms = z["__chroms__"].tolist()
            by_chrm: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
            for c in chroms:
                by_chrm[str(c)] = (
                    z[f"s::{c}"].astype(np.int64, copy=False),
                    z[f"e::{c}"].astype(np.int64, copy=False),
                    z[f"n::{c}"],
                )
        return cls(by_chrm)

    @classmethod
    def build_from_db(cls, db_path: Path) -> "_LocationIndex":
        """Scan genes.db, unpickle each BLOB, extract (chrm, gene_start, gene_end, identifier).

        Identifier prefers ``gene_name``; falls back to ``gene_id`` when the name is empty
        (common for Ensembl lncRNAs/pseudogenes).
        """
        if not db_path.exists():
            raise FileNotFoundError(f"genes.db not found at {db_path}")
        per_chrm: Dict[str, List[Tuple[int, int, str]]] = {}
        conn = sqlite3.connect(str(db_path))
        try:
            cur = conn.execute("SELECT gene_name, gene_id, data FROM genes")
            for gene_name, gene_id, blob in cur:
                if blob is None:
                    continue
                try:
                    data = pickle.loads(blob)
                except Exception:
                    continue
                chrm = data.get("chrm")
                start = data.get("gene_start")
                end = data.get("gene_end")
                if chrm is None or start is None or end is None:
                    continue
                ident = gene_name if (gene_name and str(gene_name).strip()) else gene_id
                if not ident:
                    continue
                if start > end:
                    start, end = end, start
                key = _normalize_chrm(str(chrm))
                per_chrm.setdefault(key, []).append((int(start), int(end), str(ident)))
        finally:
            conn.close()

        by_chrm: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for c, rows in per_chrm.items():
            rows.sort(key=lambda r: r[0])
            starts = np.fromiter((r[0] for r in rows), dtype=np.int64, count=len(rows))
            ends = np.fromiter((r[1] for r in rows), dtype=np.int64, count=len(rows))
            names = np.array([r[2] for r in rows], dtype=object)
            by_chrm[c] = (starts, ends, names)
        return cls(by_chrm)


def _get_index(organism: Optional[str] = None, rebuild: bool = False) -> _LocationIndex:
    """Return the cached index for organism, loading or building as needed."""
    if organism is None:
        organism = get_default_organism()
    if not rebuild and organism in _INDEX_CACHE:
        return _INDEX_CACHE[organism]

    with _BUILD_LOCK:
        if not rebuild and organism in _INDEX_CACHE:
            return _INDEX_CACHE[organism]
        npz_path = _locations_path(organism)
        if npz_path.exists() and not rebuild:
            index = _LocationIndex.load(npz_path)
        else:
            db_path = _genes_db_path(organism)
            index = _LocationIndex.build_from_db(db_path)
            try:
                index.save(npz_path)
            except OSError:
                # Read-only data dir is fine; we still have it in memory.
                pass
        _INDEX_CACHE[organism] = index
        return index


def build_location_index(organism: Optional[str] = None, force: bool = False) -> Path:
    """Build (or rebuild) the on-disk location index for an organism.

    Returns the sidecar path. Use ``force=True`` to overwrite an existing index.
    """
    if organism is None:
        organism = get_default_organism()
    npz_path = _locations_path(organism)
    if npz_path.exists() and not force:
        _INDEX_CACHE.pop(organism, None)
        _get_index(organism)
        return npz_path
    _INDEX_CACHE.pop(organism, None)
    _get_index(organism, rebuild=True)
    return npz_path


def gene_names_at_position(
    chrm: str, pos: PosArg, organism: Optional[str] = None
) -> List[str]:
    """Fast name-only lookup: gene names overlapping a point or range. No BLOB load."""
    return _get_index(organism).query_names(chrm, pos)


def write_location_index_from_rows(
    organism_dir: Path, rows: List[Tuple[str, int, int, str]]
) -> Path:
    """Write a location index directly from rows, without scanning genes.db.

    Used by the GTF build path so a fresh setup gets the sidecar for free.

    Each row is ``(chrm, gene_start, gene_end, gene_name)``.
    """
    per_chrm: Dict[str, List[Tuple[int, int, str]]] = {}
    for chrm, start, end, name in rows:
        if chrm is None or start is None or end is None or name is None:
            continue
        if start > end:
            start, end = end, start
        per_chrm.setdefault(_normalize_chrm(str(chrm)), []).append(
            (int(start), int(end), str(name))
        )
    by_chrm: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for c, rs in per_chrm.items():
        rs.sort(key=lambda r: r[0])
        starts = np.fromiter((r[0] for r in rs), dtype=np.int64, count=len(rs))
        ends = np.fromiter((r[1] for r in rs), dtype=np.int64, count=len(rs))
        names = np.array([r[2] for r in rs], dtype=object)
        by_chrm[c] = (starts, ends, names)
    index = _LocationIndex(by_chrm)
    out = Path(organism_dir) / "gene_locations.npz"
    index.save(out)
    return out
