"""Typed exceptions raised by SeqMat.

All custom exceptions inherit from :class:`SeqMatError` so callers can catch
the whole family with a single ``except SeqMatError``.
"""
from __future__ import annotations


class SeqMatError(Exception):
    """Base class for all SeqMat exceptions."""


class GeneNotFoundError(SeqMatError, LookupError):
    """Raised when a gene cannot be located in the configured organism's database."""


class OrganismNotConfiguredError(SeqMatError, ValueError):
    """Raised when an organism has no data set up (run ``seqmat setup``)."""


class DataUnavailableError(SeqMatError, RuntimeError):
    """Raised when prebuilt data (genes.db, FASTA, etc.) cannot be obtained."""


# Backward-compat alias for the previous public name.
PrebuiltDataUnavailableError = DataUnavailableError


__all__ = [
    "SeqMatError",
    "GeneNotFoundError",
    "OrganismNotConfiguredError",
    "DataUnavailableError",
    "PrebuiltDataUnavailableError",
]
