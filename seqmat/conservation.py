"""Conservation score storage: load from .pkl/.db and convert between them."""
from __future__ import annotations

import pickle
import sqlite3
from pathlib import Path
from typing import Any, Dict, Optional, Union

from ._io import unload_pickle


def load_conservation(path: Union[str, Path]) -> Dict[str, Any]:
    """Load conservation scores from a ``.db`` (SQLite) or ``.pkl`` file.

    Returns ``dict[transcript_id] -> {"scores": np.ndarray, "seq": str | bytes}``.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".db":
        conn = sqlite3.connect(str(path))
        conn.row_factory = sqlite3.Row
        try:
            cursor = conn.execute("SELECT transcript_id, data FROM conservation")
            out: Dict[str, Any] = {}
            for row in cursor:
                out[str(row["transcript_id"])] = pickle.loads(row["data"])
        finally:
            conn.close()
        return out
    data = unload_pickle(path)
    if not isinstance(data, dict):
        raise ValueError("Conservation pickle must be a dict[transcript_id] -> {scores, seq}")
    return data


def convert_conservation_pkl_to_db(
    pkl_path: Union[str, Path],
    db_path: Optional[Union[str, Path]] = None,
) -> Path:
    """Convert a conservation ``.pkl`` to ``conservation.db`` (SQLite).

    Returns the path to the written .db file.
    """
    pkl_path = Path(pkl_path)
    db_path = Path(db_path) if db_path is not None else pkl_path.with_suffix(".db")
    if db_path.exists():
        db_path.unlink()
    data = unload_pickle(pkl_path)
    if not isinstance(data, dict):
        raise ValueError("Conservation pickle must be a dict[transcript_id] -> {scores, seq}")
    conn = sqlite3.connect(str(db_path))
    try:
        conn.execute("CREATE TABLE conservation (transcript_id TEXT PRIMARY KEY, data BLOB)")
        conn.executemany(
            "INSERT INTO conservation (transcript_id, data) VALUES (?, ?)",
            ((str(tid), pickle.dumps(blob)) for tid, blob in data.items()),
        )
        conn.commit()
    finally:
        conn.close()
    return db_path


__all__ = ["load_conservation", "convert_conservation_pkl_to_db"]
