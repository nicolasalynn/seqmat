"""Transcript: RNA transcript with exons, introns, splice sites, and optional protein translation."""
from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np

_CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def _translate(seq: str) -> str:
    """Translate a nucleotide sequence to protein using the standard codon table."""
    seq = seq.upper()
    # Trim to multiple of 3 (avoids partial-codon warnings)
    trimmed = seq[:len(seq) - len(seq) % 3] if len(seq) % 3 else seq
    return ''.join(_CODON_TABLE.get(trimmed[i:i+3], 'X') for i in range(0, len(trimmed), 3))

from .seqmat import SeqMat, _fetch_fasta_region
from .config import get_organism_config, get_default_organism


class Transcript:
    """Transcript with genomic boundaries, donors/acceptors, pre_mrna, mature_mrna, and optional ORF/protein."""

    def __init__(self, d: Dict[str, Any], organism: Optional[str] = None):
        """Build Transcript from dict of attributes. Requires transcript_start, transcript_end, rev, chrm."""
        self._cons_vector = None
        array_fields = {"acceptors", "donors", "cons_vector", "rev"}
        for k, v in d.items():
            if k == "cons_vector":
                self._cons_vector = np.asarray(v, dtype=np.float32) if v is not None else None
                continue
            if k in array_fields and v is not None:
                v = np.array(v)
            setattr(self, k, v)

        self.organism = organism if organism is not None else get_default_organism()
        required_attrs = ["transcript_start", "transcript_end", "rev", "chrm"]
        missing = [attr for attr in required_attrs if not hasattr(self, attr)]
        if missing:
            raise AssertionError(f"Transcript is missing required attributes: {missing}")
        if not hasattr(self, "donors") or self.donors is None:
            self.donors = np.array([])
        if not hasattr(self, "acceptors") or self.acceptors is None:
            self.acceptors = np.array([])
        if not hasattr(self, "cons_available"):
            self.cons_available = False
        self.protein_coding = hasattr(self, "TIS") and hasattr(self, "TTS")
        self.transcript_upper = max(self.transcript_start, self.transcript_end)
        self.transcript_lower = min(self.transcript_start, self.transcript_end)
        self._pre_mrna = None
        if self.cons_available and hasattr(self, "cons_seq") and self._cons_vector is not None:
            if self.cons_seq.endswith("*") and len(self.cons_seq) == len(self._cons_vector):
                self._cons_vector = self._cons_vector[:-1]
                self.cons_seq = self.cons_seq[:-1]

    def __repr__(self) -> str:
        """Official string representation."""
        return f"Transcript({getattr(self, 'transcript_id', 'unknown_id')})"

    def __str__(self) -> str:
        """User-friendly string representation of the transcript."""
        transcript_biotype = getattr(self, 'transcript_biotype', 'unknown').replace('_', ' ').title()
        primary = getattr(self, 'primary_transcript', False)
        return f"Transcript {getattr(self, 'transcript_id', 'unknown_id')}, " \
               f"Type: {transcript_biotype}, Primary: {primary}"

    def __len__(self) -> int:
        """Length of the transcript sequence."""
        return len(getattr(self, 'transcript_seq', ''))

    def __eq__(self, other: object) -> bool:
        """Check equality of two transcripts based on their transcript sequences."""
        if not isinstance(other, Transcript):
            return NotImplemented
        return self.transcript_seq == other.transcript_seq

    def __contains__(self, subvalue: Any) -> bool:
        """Check if a given subsequence is contained within the pre_mRNA."""
        if not hasattr(subvalue, 'seq_array'):
            return False
        return np.all(np.isin(subvalue.index, self.pre_mrna.index))

    def clone(self) -> Transcript:
        """Returns a deep copy of this Transcript instance."""
        return copy.deepcopy(self)

    @property
    def pre_mrna(self) -> SeqMat:
        """Pre-mRNA sequence. Lazy-loaded from FASTA on first access."""
        if self._pre_mrna is None:
            self.generate_pre_mrna()
        return self._pre_mrna

    @pre_mrna.setter
    def pre_mrna(self, value: Optional[SeqMat]) -> None:
        self._pre_mrna = value

    @property
    def cons_vector(self) -> np.ndarray:
        """Per-residue conservation scores; ones vector when conservation data not loaded."""
        if self._cons_vector is not None:
            return self._cons_vector
        n = len(getattr(self, "protein", "")) if getattr(self, "protein", "") else len(self.pre_mrna)
        return np.ones(n, dtype=np.float32)

    @property
    def exons(self) -> List[Tuple[int, int]]:
        """Return a list of exon boundary tuples (acceptor, donor)."""
        exon_starts = np.concatenate(([self.transcript_start], self.acceptors))
        exon_ends = np.concatenate((self.donors, [self.transcript_end]))
        return list(zip(exon_starts, exon_ends))

    @property
    def exons_pos(self) -> List[Tuple[int, int]]:
        """Return exons with positions adjusted for strand orientation."""
        exon_positions = self.exons
        if self.rev:
            # Reverse order and swap coordinates for reverse strand
            exon_positions = [(end, start) for start, end in exon_positions[::-1]]
        return exon_positions

    @property
    def introns(self) -> List[Tuple[int, int]]:
        """Return a list of intron boundaries derived from donors and acceptors."""
        valid_donors = self.donors[self.donors != self.transcript_end]
        valid_acceptors = self.acceptors[self.acceptors != self.transcript_start]
        
        # Adjust intron boundaries to exclude exon splice sites
        introns = []
        for donor, acceptor in zip(valid_donors, valid_acceptors):
            if self.rev:
                # For reverse strand: intron from (donor-1) to (acceptor+1)
                intron_start = donor - 1
                intron_end = acceptor + 1
            else:
                # For forward strand: intron from (donor+1) to (acceptor-1)
                intron_start = donor + 1
                intron_end = acceptor - 1
            introns.append((intron_start, intron_end))
        
        return introns

    @property
    def introns_pos(self) -> List[Tuple[int, int]]:
        """Return introns with positions adjusted for strand orientation."""
        intron_positions = self.introns
        if self.rev:
            intron_positions = [(end, start) for start, end in intron_positions[::-1]]
        return intron_positions

    def _fix_and_check_introns(self) -> 'Transcript':
        """
        Ensure acceptors and donors are sorted and unique, and validate exon/intron structures.

        Raises:
            ValueError: If there are mismatches or ordering issues in exons/introns

        Returns:
            The current Transcript object (for chaining)
        """
        # Ensure uniqueness and correct ordering based on strand
        self.acceptors = np.unique(self.acceptors)
        self.donors = np.unique(self.donors)

        if self.rev:
            self.acceptors = np.sort(self.acceptors)[::-1]
            self.donors = np.sort(self.donors)[::-1]
        else:
            self.acceptors = np.sort(self.acceptors)
            self.donors = np.sort(self.donors)

        # Validation checks
        if self._exon_intron_matchup_flag():
            raise ValueError("Unequal number of acceptors and donors.")

        if self._exon_intron_order_flag():
            raise ValueError("Exon/intron order out of position.")

        if self._transcript_boundary_flag():
            raise ValueError("Transcript boundaries must straddle acceptors and donors.")

        return self

    def _exon_intron_matchup_flag(self) -> bool:
        """Check if acceptors and donors count match."""
        return len(self.acceptors) != len(self.donors)

    def _exon_intron_order_flag(self) -> bool:
        """Check for ordering issues in exon boundaries."""
        return any(start > end for start, end in self.exons_pos)

    def _transcript_boundary_flag(self) -> bool:
        """Check if boundaries are within the transcript start/end range."""
        if not len(self.acceptors) and not len(self.donors):
            return False
        min_boundary = np.min(np.concatenate((self.acceptors, self.donors)))
        max_boundary = np.max(np.concatenate((self.acceptors, self.donors)))
        return (self.transcript_lower > min_boundary) or (self.transcript_upper < max_boundary)

    @property
    def exonic_indices(self) -> np.ndarray:
        """Return the indices covering exons in the transcript."""
        return np.concatenate([np.arange(a, b + 1) for a, b in self.exons_pos])

    def _resolve_fasta_path(self, fasta_path: Optional[Path] = None) -> Tuple[Path, List[str]]:
        """Resolve FASTA file path and contig names for this transcript."""
        if fasta_path is None:
            config = get_organism_config(self.organism)
            fasta_path = config.get("fasta_full_genome") or config.get("fasta")
            if fasta_path is None:
                fasta_path = config.get("CHROM_SOURCE") and (Path(config["CHROM_SOURCE"]) / f"chr{self.chrm}.fasta")
            fasta_path = Path(fasta_path) if fasta_path else None
            if not fasta_path or not fasta_path.exists():
                raise ValueError(
                    "No FASTA path in config (fasta_full_genome or fasta). "
                    "Run setup_genomics_data() to install the full-genome FASTA (e.g. hg38.fa)."
                )
        raw_chr = str(self.chrm).replace("chr", "") or "17"
        with_chr = f"chr{raw_chr}"
        contigs = [with_chr, raw_chr] if with_chr != raw_chr else [with_chr]
        return Path(fasta_path), contigs

    def _assemble_mature_from_fasta(self) -> SeqMat:
        """Fetch exons directly from FASTA and assemble mature mRNA (skips pre-mRNA)."""
        fasta_path, contigs = self._resolve_fasta_path()

        # Sort exons by ascending genomic position
        exon_intervals = sorted((min(a, b), max(a, b)) for a, b in self.exons)

        last_err = None
        for contig in contigs:
            try:
                seqs = []
                all_indices = []
                for lo, hi in exon_intervals:
                    seq_str = _fetch_fasta_region(str(fasta_path), contig, lo, hi)
                    seqs.append(seq_str)
                    all_indices.append(np.arange(lo, hi + 1, dtype=np.int64))

                mature = SeqMat(''.join(seqs), indices=np.concatenate(all_indices))
                if self.rev:
                    mature.reverse_complement()
                return mature
            except Exception as e:
                last_err = e
        if last_err is not None:
            raise last_err
        raise ValueError(f"Could not fetch exons from FASTA for {self.chrm}")

    def pull_pre_mrna_from_fasta(
        self,
        fasta_path: Optional[Path] = None,
        upstream: int = 0,
        downstream: int = 0,
        region_start: Optional[int] = None,
        region_end: Optional[int] = None,
    ) -> Dict[str, Any]:
        """Pre-mRNA sequence from the single full-genome FASTA (e.g. hg38.fa). Optional upstream/downstream."""
        fasta_path, contigs_to_try = self._resolve_fasta_path(fasta_path)
        if region_start is not None or region_end is not None:
            start = region_start if region_start is not None else max(1, self.transcript_lower - upstream)
            end = region_end if region_end is not None else (self.transcript_upper + downstream)
        else:
            start = max(1, self.transcript_lower - upstream)
            end = self.transcript_upper + downstream
        last_err = None
        for contig in contigs_to_try:
            try:
                seq_mat = SeqMat.from_fasta_file(fasta_path, contig, start, end)
                if len(seq_mat.seq) > 0:
                    return {"seq": seq_mat.seq, "indices": seq_mat.index}
            except Exception as e:
                last_err = e
        if last_err is not None:
            raise last_err
        raise ValueError(f"Could not fetch {self.chrm}:{start}-{end} from FASTA (tried contigs {contigs_to_try})")

    def generate_pre_mrna(
        self,
        upstream: int = 0,
        downstream: int = 0,
        region_start: Optional[int] = None,
        region_end: Optional[int] = None,
    ) -> "Transcript":
        """Load pre-mRNA and set self.pre_mrna. Use upstream/downstream to add flanking (bp).

        Optional region_start/region_end restrict the FASTA read to a genomic
        sub-interval (e.g. ±7500 bp around a mutation site) instead of fetching
        the entire transcript, which can be orders of magnitude faster for large
        genes.
        """
        try:
            seq_data = self.pull_pre_mrna_from_fasta(
                upstream=upstream, downstream=downstream,
                region_start=region_start, region_end=region_end,
            )
            pre_mrna = SeqMat(**seq_data)
        except Exception:
            if region_start is not None or region_end is not None:
                lo = region_start if region_start is not None else max(1, self.transcript_lower - upstream)
                hi = region_end if region_end is not None else (self.transcript_upper + downstream)
            else:
                lo = max(1, self.transcript_lower - upstream)
                hi = self.transcript_upper + downstream
            length = hi - lo + 1
            indices = np.arange(lo, lo + length)
            pre_mrna = SeqMat("N" * length, indices=indices)
        if self.rev:
            pre_mrna.reverse_complement()
        self.pre_mrna = pre_mrna
        return self

    def generate_mature_mrna(self, inplace: bool = True) -> Union['Transcript', SeqMat]:
        """
        Generate the mature mRNA by concatenating exon regions from pre_mRNA.

        When pre-mRNA has not been loaded yet, exons are fetched directly from
        the indexed FASTA — skipping the full-gene read entirely.

        Args:
            inplace: If True, set self.mature_mrna, else return a new SeqMat

        Returns:
            The Transcript object (if inplace=True) or a SeqMat (if inplace=False)
        """
        self._fix_and_check_introns()

        if self._pre_mrna is not None:
            # Pre-mRNA loaded (possibly mutated): splice from it
            result = self.pre_mrna.remove_regions(self.introns)
        else:
            # Fast path: fetch exons directly, skip full-gene read
            try:
                result = self._assemble_mature_from_fasta()
            except Exception:
                # Fall back to full pre-mRNA load + splice
                result = self.pre_mrna.remove_regions(self.introns)

        if inplace:
            self.mature_mrna = result
            return self
        return result

    @property
    def orf(self) -> Union[SeqMat, 'Transcript']:
        """
        Return the ORF (Open Reading Frame) SeqMat object, if TIS and TTS are available.

        Returns:
            The ORF SeqMat if TIS/TTS are set, else self
        """
        if not self.protein_coding:
            print("Cannot create protein without set TIS and TTS values.")
            return self

        # Extract ORF region from mature mRNA via clone (no string re-parse)
        if hasattr(self, 'mature_mrna') and self.mature_mrna is not None:
            return self.mature_mrna.clone(self.TIS, self.TTS)

        return self

    def generate_protein(self, inplace: bool = True) -> Union['Transcript', str]:
        """
        Translate the ORF into a protein sequence.

        Args:
            inplace: If True, store protein in self. Otherwise, return it

        Returns:
            The Transcript object if inplace=True, else the protein sequence
        """
        if not self.protein_coding:
            return self if inplace else ""

        # Translate the ORF to protein
        orf_seq = self.orf
        if isinstance(orf_seq, SeqMat):
            protein = _translate(orf_seq.seq).strip('*')
        else:
            protein = ""

        if inplace:
            self.protein = protein
            if self._cons_vector is not None and len(self._cons_vector) != len(protein):
                self._cons_vector = np.ones(len(protein), dtype=np.float32)
            return self
        
        return protein