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

def get_default_organism() -> str:
    """Get the default organism from config or fallback"""
    config = load_config()
    return config.get('default_organism', DEFAULT_SETTINGS['default_organism'])

def get_available_organisms() -> List[str]:
    """Get list of available organisms from config and defaults"""
    config = load_config()
    configured_organisms = set(config.keys()) - {'default_organism', 'directory_structure'}
    default_organisms = set(DEFAULT_ORGANISM_DATA.keys())
    return sorted(configured_organisms | default_organisms)

def get_organism_info(organism: str) -> Dict[str, Any]:
    """Get organism information including name and URLs"""
    config = load_config()
    
    if organism in config and isinstance(config[organism], dict):
        org_config = config[organism]
        # Merge with defaults if available
        if organism in DEFAULT_ORGANISM_DATA:
            default_data = DEFAULT_ORGANISM_DATA[organism].copy()
            default_data.update(org_config)
            return default_data
        return org_config
    elif organism in DEFAULT_ORGANISM_DATA:
        return DEFAULT_ORGANISM_DATA[organism]
    else:
        raise ValueError(f"Organism '{organism}' not configured. Available: {get_available_organisms()}")

def get_organism_config(organism: Optional[str] = None) -> Dict[str, Path]:
    """Get configuration paths for a specific organism"""
    if organism is None:
        organism = get_default_organism()
        
    config = load_config()
    if organism not in config:
        raise ValueError(f"Organism '{organism}' not configured. Run setup_genomics_data() first.")
    
    # Convert string paths to Path objects
    org_config = config[organism]
    
    # Handle case where org_config might be a string instead of dict
    if isinstance(org_config, str):
        raise ValueError(f"Invalid configuration for organism '{organism}'. "
                        f"Expected dictionary but got string: {org_config}")
    
    if not isinstance(org_config, dict):
        raise ValueError(f"Invalid configuration for organism '{organism}'. "
                        f"Expected dictionary but got {type(org_config)}")
    
    return {k: Path(v) for k, v in org_config.items() if isinstance(v, str)}

def get_directory_config() -> Dict[str, str]:
    """Get directory structure configuration"""
    config = load_config()
    return config.get('directory_structure', DEFAULT_SETTINGS['directory_structure'])