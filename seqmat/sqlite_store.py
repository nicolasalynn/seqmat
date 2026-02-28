"""SQLite backend for gene data loading.

Setup writes genes.db in the organism directory (or it is downloaded from S3).
Gene.from_file() tries SQLite first (after LMDB), then falls back to per-gene pickle files.
Single indexed DB for fast lookups.
"""

import pickle
from pathlib import Path
from typing import Any, Dict, Optional

from .config import get_organism_config, get_default_organism, load_config


def get_genes_db_path(organism: Optional[str] = None) -> Optional[str]:
    """Return genes_db path from config if set and file exists; else None."""
    if organism is None:
        organism = get_default_organism()
    try:
        config = get_organism_config(organism)
    except ValueError:
        return None
    path = config.get("genes_db")
    if not path:
        return None
    p = Path(path)
    return str(p) if p.exists() else None


def load_gene_from_sqlite(gene_name: str, organism: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Load gene dict from genes.db if configured; otherwise return None."""
    db_path = get_genes_db_path(organism)
    if not db_path:
        return None
    try:
        import sqlite3
        conn = sqlite3.connect(db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.execute(
            "SELECT data FROM genes WHERE UPPER(gene_name) = UPPER(?) LIMIT 1",
            (gene_name.strip(),),
        )
        row = cursor.fetchone()
        conn.close()
        if row is None:
            return None
        return pickle.loads(row["data"])
    except Exception:
        return None
