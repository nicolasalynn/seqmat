"""
SeqMat - Lightning-fast genomic sequence matrix library

A comprehensive Python library for genomic sequence analysis with full mutation tracking,
splicing analysis, and sequence manipulation.
"""

__version__ = "0.1.3"
__author__ = "Nicolas Lynn Vila"
__email__ = "nicolasalynn@gmail.com"

from .seqmat import SeqMat
from .gene import Gene
from .transcript import Transcript
from .utils import (
    setup_genomics_data, 
    load_config, 
    save_config,
    list_available_organisms,
    list_supported_organisms,
    get_organism_info,
    list_gene_biotypes,
    count_genes,
    get_gene_list,
    data_summary,
    print_data_summary,
    search_genes
)

__all__ = [
    "SeqMat",
    "Gene", 
    "Transcript",
    "setup_genomics_data",
    "load_config",
    "save_config",
    "list_available_organisms",
    "list_supported_organisms", 
    "get_organism_info",
    "list_gene_biotypes",
    "count_genes",
    "get_gene_list",
    "data_summary",
    "print_data_summary",
    "search_genes"
]