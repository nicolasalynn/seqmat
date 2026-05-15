"""SeqMat — fast, vectorized genomic sequences with first-class mutation tracking."""

__version__ = "1.5.0"
__author__ = "Nicolas Lynn Vila"
__email__ = "nicolasalynn@gmail.com"

# Core classes
from .seqmat import SeqMat
from .gene import Gene
from .transcript import Transcript

# Config / paths
from .config import (
    get_config_dir,
    get_config_file,
    get_data_base,
    get_data_dir,
    get_default_organism,
    load_config,
    save_config,
)

# Optional LMDB backend
from .lmdb_store import build_lmdb

# Position → gene lookup
from .locator import build_location_index, gene_names_at_position

# Data setup
from .data_setup import set_fasta_path, setup_genomics_data

# Discovery
from .discovery import (
    available_genes,
    count_genes,
    data_summary,
    get_all_genes,
    get_gene_list,
    get_organism_info,
    list_available_organisms,
    list_gene_biotypes,
    list_supported_organisms,
    print_data_summary,
    search_genes,
)

# Errors
from .errors import (
    DataUnavailableError,
    GeneNotFoundError,
    OrganismNotConfiguredError,
    PrebuiltDataUnavailableError,
    SeqMatError,
)

__all__ = [
    # Core
    "SeqMat", "Gene", "Transcript",
    # Config / paths
    "get_default_organism", "get_data_dir", "get_config_dir",
    "get_config_file", "get_data_base", "load_config", "save_config",
    # Setup
    "setup_genomics_data", "set_fasta_path",
    # Discovery
    "list_available_organisms", "list_supported_organisms", "get_organism_info",
    "list_gene_biotypes", "count_genes", "get_gene_list", "get_all_genes",
    "data_summary", "print_data_summary", "search_genes", "available_genes",
    # Position lookup
    "gene_names_at_position", "build_location_index",
    # LMDB
    "build_lmdb",
    # Errors
    "SeqMatError", "GeneNotFoundError", "OrganismNotConfiguredError",
    "DataUnavailableError", "PrebuiltDataUnavailableError",
]
