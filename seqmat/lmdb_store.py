"""LMDB backend for fast gene data loading.

Optional: when the ``lmdb`` package is installed and ``gene_lmdb_path`` is set,
Gene.from_file() will try LMDB first and fall back to per-gene pickle files otherwise.
All behavior is unchanged when LMDB is not used.
"""

import pickle
import shutil
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from .config import get_organism_config, get_default_organism, load_config

try:
    import lmdb as _lmdb
except ImportError:
    _lmdb = None  # type: ignore[assignment]

_ENV_CACHE: Dict[str, Any] = {}


def _require_lmdb() -> None:
    """Raise ImportError if lmdb is not installed."""
    if _lmdb is None:
        raise ImportError(
            "The 'lmdb' package is required for LMDB support. "
            "Install with: pip install seqmat[lmdb]"
        )


def get_lmdb_config(organism: Optional[str] = None) -> Tuple[Optional[str], bool]:
    """Return (lmdb_path, local_staging) from config. (None, False) when not set."""
    config = load_config()
    if organism is None:
        organism = get_default_organism()
    org_conf = config.get(organism, {})
    if isinstance(org_conf, dict):
        lmdb_path = org_conf.get("gene_lmdb_path")
        if lmdb_path is not None:
            staging = org_conf.get("gene_lmdb_local_staging", False)
            return str(lmdb_path), bool(staging)
    lmdb_path = config.get("gene_lmdb_path")
    staging = config.get("gene_lmdb_local_staging", False)
    return (str(lmdb_path) if lmdb_path else None), bool(staging)


def _open_lmdb(path: str, local_staging: bool = False) -> Any:
    """Open read-only LMDB at path; optionally copy to /tmp first. Cached per path."""
    _require_lmdb()
    if local_staging:
        import tempfile
        staging_dir = Path(tempfile.gettempdir()) / f"seqmat_lmdb_{Path(path).name}"
        if not staging_dir.exists():
            shutil.copytree(path, str(staging_dir))
        path = str(staging_dir)
    if path in _ENV_CACHE:
        return _ENV_CACHE[path]
    env = _lmdb.open(
        path,
        readonly=True,
        readahead=True,
        lock=False,
        max_dbs=0,
        map_size=0,
    )
    _ENV_CACHE[path] = env
    return env


def load_gene_from_lmdb(gene_name: str, organism: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """Load gene dict from LMDB if configured; otherwise return None. No side effects."""
    if _lmdb is None:
        return None
    if organism is None:
        organism = get_default_organism()
    lmdb_path, local_staging = get_lmdb_config(organism)
    if lmdb_path is None:
        return None
    try:
        env = _open_lmdb(lmdb_path, local_staging=local_staging)
    except Exception:
        return None
    with env.begin() as txn:
        raw = txn.get(gene_name.encode("utf-8"))
        if raw is None:
            return None
        return pickle.loads(raw)


def _gene_name_from_pkl_stem(stem: str) -> str:
    """Derive LMDB key from pickle filename stem. Matches glob *_{gene_name}.pkl behavior."""
    parts = stem.split("_")
    if len(parts) >= 3 and parts[0] == "mrnas":
        return "_".join(parts[2:])
    if len(parts) >= 2:
        return "_".join(parts[1:])
    return stem


def build_lmdb(
    annotations_dir: Optional[str] = None,
    output_path: Optional[str] = None,
    organism: Optional[str] = None,
) -> str:
    """Build LMDB from per-gene pickle files. Requires pip install seqmat[lmdb]."""
    _require_lmdb()
    if organism is None:
        organism = get_default_organism()
    if annotations_dir is None:
        config = get_organism_config(organism)
        annotations_dir = str(config["MRNA_PATH"])
    ann_path = Path(annotations_dir)
    if not ann_path.exists():
        raise FileNotFoundError(f"Annotations directory not found: {ann_path}")
    if output_path is None:
        output_path = str(ann_path / "genes.lmdb")
    out = Path(output_path)
    if out.exists():
        shutil.rmtree(str(out))
    pkl_files = sorted(ann_path.glob("**/*.pkl"))
    if not pkl_files:
        raise FileNotFoundError(f"No .pkl files found under {ann_path}")
    total_bytes = sum(f.stat().st_size for f in pkl_files)
    map_size = int(total_bytes * 1.5) + 10 * 1024 * 1024
    env = _lmdb.open(str(out), map_size=map_size, max_dbs=0)
    genes_written = 0
    skipped = 0
    total_size = 0
    with env.begin(write=True) as txn:
        for pkl_file in pkl_files:
            try:
                raw_bytes = pkl_file.read_bytes()
                gene_name = _gene_name_from_pkl_stem(pkl_file.stem)
                txn.put(gene_name.encode("utf-8"), raw_bytes)
                genes_written += 1
                total_size += len(raw_bytes)
            except Exception as exc:
                print(f"  Skipped {pkl_file.name}: {exc}")
                skipped += 1
    env.close()
    print(f"LMDB built at: {out}")
    print(f"  Genes written: {genes_written:,}")
    print(f"  Total size:    {total_size / (1024 * 1024):.1f} MB")
    if skipped:
        print(f"  Skipped:       {skipped}")
    return str(out)
