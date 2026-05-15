"""Backward-compat re-export shim.

The previous monolithic ``utils.py`` has been split into focused modules:

- :mod:`seqmat._io` — pickle/JSON I/O
- :mod:`seqmat.download` — HTTP downloads and prebuilt fetch
- :mod:`seqmat.gtf` — GTF parsing and ``genes.db`` build pipeline
- :mod:`seqmat.conservation` — conservation file conversion
- :mod:`seqmat.data_setup` — :func:`setup_genomics_data`, :func:`set_fasta_path`
- :mod:`seqmat.discovery` — list/count/search/summary helpers

Everything that previously lived in :mod:`seqmat.utils` is still importable from
here (``from seqmat.utils import X``). New code should import from the specific
module instead.
"""
from __future__ import annotations

import logging
import random

# Public re-exports ----------------------------------------------------------

from ._io import dump_json, dump_pickle, unload_json, unload_pickle  # noqa: F401
from .conservation import (  # noqa: F401
    convert_conservation_pkl_to_db,
    load_conservation,
)
from .data_setup import set_fasta_path, setup_genomics_data  # noqa: F401
from .discovery import (  # noqa: F401
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
from .download import (  # noqa: F401
    PrebuiltDataUnavailableError,
    download,
    download_and_ungzip,
    download_genome_data,
    download_prebuilt_data,
)
from .errors import (  # noqa: F401
    DataUnavailableError,
    GeneNotFoundError,
    OrganismNotConfiguredError,
    SeqMatError,
)
from .gtf import (  # noqa: F401
    process_transcript,
    retrieve_and_parse_ensembl_annotations,
    split_fasta,
)

_log = logging.getLogger(__name__)


def generate_random_sequence(k: int = 15_000) -> str:
    """Return a random DNA sequence of length ``k`` over ``ACGT``. Useful for tests."""
    return "".join(random.choices("ACGT", k=k))


def test_installation(organism: str = None, verbose: bool = True) -> "dict":
    """Smoke-test a SeqMat installation against a configured organism.

    Returns a dict with ``tests_passed``, ``tests_failed``, ``errors``, ``warnings``.
    Hits the full stack: config, data paths, gene listing, ``Gene.from_file``,
    ``Gene.from_position``, FASTA read, and ``SeqMat`` mutation/complement.
    """
    from pathlib import Path

    from .config import get_default_organism, get_organism_config
    from .gene import Gene
    from .seqmat import SeqMat

    if organism is None:
        organism = get_default_organism()

    results = {"organism": organism, "tests_passed": 0, "tests_failed": 0, "errors": [], "warnings": []}

    def step(name: str, fn):
        if verbose:
            print(f"  {name}...", end=" ", flush=True)
        try:
            value = fn()
            results["tests_passed"] += 1
            if verbose:
                print("ok")
            return True, value
        except Exception as e:
            results["tests_failed"] += 1
            results["errors"].append(f"{name}: {e}")
            if verbose:
                print(f"FAIL ({e})")
            return False, None

    if verbose:
        print(f"Running comprehensive tests for {organism}")
        print("-" * 50)

    ok, config = step("Configuration", lambda: get_organism_config(organism))
    if not ok:
        return results

    step("Genes table", lambda: count_genes(organism))
    step("Gene search", lambda: search_genes(organism, "KRAS", limit=1) or (_ for _ in ()).throw(AssertionError("no match")))

    def _gene_loading():
        names = get_gene_list(organism, "protein_coding", limit=5)
        if not names:
            raise ValueError("no protein-coding genes")
        gene = Gene.from_file(names[0], organism=organism)
        if gene is None:
            raise ValueError(f"failed to load {names[0]}")
        return gene

    ok, gene = step("Gene loading", _gene_loading)
    if ok and gene is not None:
        def _position():
            hits = Gene.from_position(gene.chrm, (1, 10**9), organism=organism)
            if not hits:
                raise ValueError("Gene.from_position returned no hits")
            return hits

        step("Position lookup", _position)

    def _seqmat():
        s = SeqMat("ATCGATCGATCG")
        assert s.complement().seq == "TAGCTAGCTAGC"
        s.apply_mutations([(3, "C", "G")])
        assert s.mutations
        return s

    step("Sequence operations", _seqmat)

    if verbose:
        print("-" * 50)
        print(f"Passed: {results['tests_passed']}, failed: {results['tests_failed']}")
    return results


__all__ = [
    # I/O
    "dump_pickle", "unload_pickle", "dump_json", "unload_json",
    # Download
    "download", "download_and_ungzip",
    "download_genome_data", "download_prebuilt_data",
    "PrebuiltDataUnavailableError",
    # GTF / build
    "retrieve_and_parse_ensembl_annotations", "process_transcript", "split_fasta",
    # Conservation
    "load_conservation", "convert_conservation_pkl_to_db",
    # Setup
    "setup_genomics_data", "set_fasta_path",
    # Discovery
    "list_available_organisms", "list_supported_organisms", "get_organism_info",
    "list_gene_biotypes", "count_genes", "get_gene_list", "get_all_genes",
    "data_summary", "print_data_summary", "search_genes", "available_genes",
    # Errors
    "SeqMatError", "GeneNotFoundError", "OrganismNotConfiguredError", "DataUnavailableError",
    # Misc
    "generate_random_sequence", "test_installation",
]
