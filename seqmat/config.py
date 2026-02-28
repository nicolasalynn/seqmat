"""Configuration management for SeqMat.

Config-less mode (no config file):
  - SEQMAT_DATA_DIR  → root path where organism dirs live (e.g. /path/to/seqmat_base).
    All paths are derived: {SEQMAT_DATA_DIR}/{organism}/genes.db, *.fa, etc.
    Set SEQMAT_DEFAULT_ORGANISM for default (default: hg38). No config file needed.

Config file mode (when SEQMAT_DATA_DIR is not set):
  1. SEQMAT_CONFIG_FILE  → exact path to config.json
  2. SEQMAT_CONFIG_DIR   → directory; config is {SEQMAT_CONFIG_DIR}/config.json
  3. Default             → platformdirs user config dir
"""
import json
import os
from pathlib import Path
from typing import Dict, Any, List, Optional

from platformdirs import user_config_dir, user_data_dir


def _resolve_config_file() -> Path:
    """Single place for config file path. Env overrides, else default dir + config.json."""
    env_file = os.environ.get("SEQMAT_CONFIG_FILE", "").strip()
    if env_file:
        return Path(env_file)
    env_dir = os.environ.get("SEQMAT_CONFIG_DIR", "").strip()
    if env_dir:
        return Path(env_dir) / "config.json"
    return Path(user_config_dir("seqmat", appauthor=False)) / "config.json"


CONFIG_FILE = _resolve_config_file()
DEFAULT_CONFIG_DIR = CONFIG_FILE.parent
DEFAULT_DATA_DIR = Path(user_data_dir("seqmat", appauthor=False))

# S3 bucket for prebuilt data (genes.db, FASTA). Used by default for pip installs (no rebuild).
# Layout: {base}/{organism}/genes.db, {organism}/{organism}.fa.gz, etc.
# Bucket region is eu-north-1; use that endpoint to avoid 301 redirect (which returns XML body, not data).
# Override via env SEQMAT_PREBUILT_DATA_BASE_URL or config key "prebuilt_data_base_url".
_PREBUILT_DATA_BASE_URL_DEFAULT = "https://seqmat-prebuilt-public.s3.eu-north-1.amazonaws.com"


def get_prebuilt_data_base_url() -> str:
    """
    Base URL for prebuilt data (genes.db, FASTA, conservation). Used for both
    default pip-install downloads and build-from-sources conservation.
    Resolution: SEQMAT_PREBUILT_DATA_BASE_URL env > config prebuilt_data_base_url > default (new bucket).
    """
    env = (os.environ.get("SEQMAT_PREBUILT_DATA_BASE_URL") or "").strip()
    if env:
        return env.rstrip("/")
    config = load_config()
    url = (config.get("prebuilt_data_base_url") or "").strip()
    if url:
        return url.rstrip("/")
    return _PREBUILT_DATA_BASE_URL_DEFAULT.rstrip("/")

# Conservation for build-from-sources: use prebuilt bucket (same as default download).
# get_prebuilt_data_base_url() + f"/{organism}/conservation.pkl" or .db in utils.
DEFAULT_ORGANISM_DATA = {
    'hg38': {
        'name': 'Homo sapiens (Human)',
        'urls': {
            'fasta': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz',
            'gtf': 'https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz',
            'gtex': 'https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz'
        }
    },
    'mm39': {
        'name': 'Mus musculus (Mouse)',
        'urls': {
            'fasta': 'https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz',
            'gtf': 'https://ftp.ensembl.org/pub/release-112/gtf/mus_musculus/Mus_musculus.GRCm39.112.gtf.gz'
        }
    }
}

DEFAULT_SETTINGS = {
    'default_organism': 'hg38',
    'directory_structure': {
        'chromosomes': 'chromosomes',
        'annotations': 'annotations'
    },
    'gene_lmdb_path': None,
    'gene_lmdb_local_staging': False,
    'prebuilt_data_base_url': None,
}

def load_config() -> Dict[str, Any]:
    """Load config from CONFIG_FILE, merged with DEFAULT_SETTINGS."""
    if CONFIG_FILE.exists():
        with open(CONFIG_FILE, "r") as f:
            config = json.load(f)
            merged_config = DEFAULT_SETTINGS.copy()
            merged_config.update(config)
            return merged_config
    return DEFAULT_SETTINGS.copy()

def save_config(config: Dict[str, Any]) -> None:
    """Save config to CONFIG_FILE (creates parent dir if needed)."""
    CONFIG_FILE.parent.mkdir(parents=True, exist_ok=True)
    with open(CONFIG_FILE, "w") as f:
        json.dump(config, f, indent=2)

def get_data_base() -> Optional[Path]:
    """If SEQMAT_DATA_DIR is set and exists, return that path (config-less mode). Else None."""
    env = os.environ.get("SEQMAT_DATA_DIR", "").strip()
    if not env:
        return None
    p = Path(env).expanduser().resolve()
    return p if p.exists() else None


def get_default_organism() -> str:
    """Default organism: from SEQMAT_DEFAULT_ORGANISM in config-less mode, else from config file."""
    if get_data_base() is not None:
        return os.environ.get("SEQMAT_DEFAULT_ORGANISM", "").strip() or "hg38"
    config = load_config()
    return config.get('default_organism', DEFAULT_SETTINGS['default_organism'])

NON_ORGANISM_KEYS = {"default_organism", "directory_structure", "gene_lmdb_path", "gene_lmdb_local_staging", "prebuilt_data_base_url"}


def get_available_organisms() -> List[str]:
    """Organism IDs: in config-less mode, subdirs of SEQMAT_DATA_DIR that contain genes.db; else from config."""
    base = get_data_base()
    if base is not None:
        organisms = []
        for d in base.iterdir():
            if d.is_dir() and (d / "genes.db").exists():
                organisms.append(d.name)
        return sorted(organisms) if organisms else list(DEFAULT_ORGANISM_DATA.keys())
    config = load_config()
    configured_organisms = set(config.keys()) - NON_ORGANISM_KEYS
    default_organisms = set(DEFAULT_ORGANISM_DATA.keys())
    return sorted(configured_organisms | default_organisms)


def get_organism_info(organism: str) -> Dict[str, Any]:
    """Organism metadata (name, URLs) from DEFAULT_ORGANISM_DATA; paths from get_organism_config in config-less mode."""
    if organism in DEFAULT_ORGANISM_DATA:
        return DEFAULT_ORGANISM_DATA[organism].copy()
    if get_data_base() is not None and organism in get_available_organisms():
        return {"name": organism, "urls": {}}
    config = load_config()
    if organism in config and isinstance(config[organism], dict):
        org_config = config[organism]
        if organism in DEFAULT_ORGANISM_DATA:
            out = DEFAULT_ORGANISM_DATA[organism].copy()
            out.update(org_config)
            return out
        return org_config
    raise ValueError(f"Organism '{organism}' not configured. Available: {get_available_organisms()}")


def _organism_config_from_data_base(organism: str, base: Path) -> Dict[str, Path]:
    """Build path dict from data base + organism (no config file). Hardcoded layout."""
    org_path = base / organism
    if not org_path.exists():
        raise ValueError(f"Organism '{organism}' not found under {base}. Run setup_genomics_data() first.")
    # Fasta: single .fa file in organism dir (from setup)
    fa_files = list(org_path.glob("*.fa"))
    config = {
        "BASE": org_path,
        "MRNA_PATH": org_path,
        "CHROM_SOURCE": org_path,
        "genes_db": org_path / "genes.db",
        "MISSPLICING_PATH": org_path / "missplicing",
        "ONCOSPLICE_PATH": org_path / "oncosplice",
        "TEMP": org_path / "temp",
    }
    if fa_files:
        config["fasta_full_genome"] = fa_files[0]
    return config


def get_organism_config(organism: Optional[str] = None) -> Dict[str, Path]:
    """Return path config for an organism. In config-less mode (SEQMAT_DATA_DIR), paths are derived from it."""
    if organism is None:
        organism = get_default_organism()
    base = get_data_base()
    if base is not None:
        return _organism_config_from_data_base(organism, base)
    config = load_config()
    if organism not in config:
        raise ValueError(f"Organism '{organism}' not configured. Run setup_genomics_data() first.")
    org_config = config[organism]
    if isinstance(org_config, str):
        raise ValueError(f"Invalid configuration for organism '{organism}'. Expected dictionary.")
    if not isinstance(org_config, dict):
        raise ValueError(f"Invalid configuration for organism '{organism}'.")
    out = {k: Path(v) for k, v in org_config.items() if isinstance(v, str)}
    # Prefer genes.db when it exists: if not in config, infer from BASE or MRNA_PATH (same dir as setup)
    if "genes_db" not in out:
        base_or_mrna = out.get("BASE") or out.get("MRNA_PATH")
        if base_or_mrna is not None:
            out["genes_db"] = Path(base_or_mrna) / "genes.db"
    return out

def get_directory_config() -> Dict[str, str]:
    """Return directory_structure (e.g. chromosomes, annotations folder names)."""
    config = load_config()
    return config.get('directory_structure', DEFAULT_SETTINGS['directory_structure'])

def get_data_dir() -> Path:
    """Return user data directory for seqmat (platformdirs)."""
    return DEFAULT_DATA_DIR


def get_config_dir() -> Path:
    """Return the directory containing the active config file."""
    return DEFAULT_CONFIG_DIR


def get_config_file() -> Path:
    """Return the active config file path (for display or override)."""
    return CONFIG_FILE