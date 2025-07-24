"""Configuration management for SeqMat"""
import os
import json
from pathlib import Path
from typing import Dict, Any, List, Optional

DEFAULT_CONFIG_DIR = Path.home() / '.seqmat'
CONFIG_FILE = DEFAULT_CONFIG_DIR / 'config.json'

# Default organism data sources - can be overridden in config
DEFAULT_ORGANISM_DATA = {
    'hg38': {
        'name': 'Homo sapiens (Human)',
        'urls': {
            'fasta': 'https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz',
            'gtf': 'https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz',
            'conservation': 'https://genome-data-public-access.s3.eu-north-1.amazonaws.com/conservation.pkl',
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
    }
}

def load_config() -> Dict[str, Any]:
    """Load configuration from user's home directory"""
    if CONFIG_FILE.exists():
        with open(CONFIG_FILE, 'r') as f:
            config = json.load(f)
            # Merge with default settings
            merged_config = DEFAULT_SETTINGS.copy()
            merged_config.update(config)
            return merged_config
    return DEFAULT_SETTINGS.copy()

def save_config(config: Dict[str, Any]) -> None:
    """Save configuration to user's home directory"""
    DEFAULT_CONFIG_DIR.mkdir(exist_ok=True)
    with open(CONFIG_FILE, 'w') as f:
        json.dump(config, f, indent=2)

def get_organism_config(organism: str = 'hg38') -> Dict[str, Path]:
    """Get configuration paths for a specific organism"""
    config = load_config()
    if organism not in config:
        raise ValueError(f"Organism '{organism}' not configured. Run setup_genomics_data() first.")
    
    # Convert string paths to Path objects
    org_config = config[organism]
    return {k: Path(v) for k, v in org_config.items()}