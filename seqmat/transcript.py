"""Transcript class for representing RNA transcripts with genomic information"""
from __future__ import annotations
from typing import Any, Optional, Union, Dict, List, Tuple
import numpy as np
import copy
from pathlib import Path

try:
    from Bio.Seq import Seq
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

from .seqmat import SeqMat
from .config import get_organism_config, get_default_organism
from .utils import unload_pickle


class Transcript:
    """
    Represents a transcript with associated genomic information such as exons, introns, and sequences.

    Attributes include:
    - transcript_start, transcript_end: Transcript boundaries
    - rev: Whether transcript is on reverse strand
    - chrm: Chromosome
    - donors, acceptors: Splice sites
    - cons_vector, cons_seq: Conservation data
    - transcript_seq: Transcript sequence
    - transcript_biotype: Type of transcript
    - primary_transcript: Whether this is the primary transcript
    - transcript_id: Unique identifier
    - TIS, TTS: Translation initiation/termination sites (if protein-coding)
    """

    def __init__(self, d: Dict[str, Any], organism: Optional[str] = None):
        """
        Initialize a Transcript object from a dictionary of attributes.

        Args:
            d: Dictionary containing transcript attributes and data
            organism: Genome build or organism reference (e.g., 'hg38')

        Raises:
            AssertionError: If required attributes are missing
        """
        # Convert certain attributes to NumPy arrays for consistent processing
        array_fields = {'acceptors', 'donors', 'cons_vector', 'rev'}
        for k, v in d.items():
            if k in array_fields and v is not None:
                v = np.array(v)
            setattr(self, k, v)

        self.organism: str = organism if organism is not None else get_default_organism()

        # Required attributes to form a valid transcript object
        required_attrs = ['transcript_start', 'transcript_end', 'rev', 'chrm']
        missing = [attr for attr in required_attrs if not hasattr(self, attr)]
        if missing:
            raise AssertionError(f"Transcript is missing required attributes: {missing}")

        # Default fallback values for optional attributes
        if not hasattr(self, 'donors') or self.donors is None:
            self.donors = np.array([])
        if not hasattr(self, 'acceptors') or self.acceptors is None:
            self.acceptors = np.array([])
        if not hasattr(self, 'cons_available'):
            self.cons_available = False

        # Determine if transcript is protein-coding
        self.protein_coding: bool = hasattr(self, 'TIS') and hasattr(self, 'TTS')

        # Calculate transcript boundaries
        self.transcript_upper = max(self.transcript_start, self.transcript_end)
        self.transcript_lower = min(self.transcript_start, self.transcript_end)

        # Generate pre-mRNA sequence data
        self.generate_pre_mrna()

        # If consensus data is available and ends with '*', adjust cons_vector and cons_seq
        if self.cons_available and hasattr(self, 'cons_seq') and hasattr(self, 'cons_vector'):
            if self.cons_seq.endswith('*') and len(self.cons_seq) == len(self.cons_vector):
                self.cons_vector = self.cons_vector[:-1]
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

    def pull_pre_mrna_from_fasta(self, fasta_path: Optional[Path] = None) -> Dict[str, Any]:
        """
        Retrieve the pre-mRNA sequence from a FASTA file.
        
        Args:
            fasta_path: Optional path to FASTA file. If None, uses config path.
            
        Returns:
            Dictionary with 'seq' and 'indices' keys
        """
        if fasta_path is None:
            config = get_organism_config(self.organism)
            fasta_path = config['CHROM_SOURCE'] / f'chr{self.chrm}.fasta'
        
        # Read the sequence directly
        seq_mat = SeqMat.from_fasta_file(
            fasta_path,
            f'chr{self.chrm}',
            self.transcript_lower,
            self.transcript_upper
        )
        
        return {
            'seq': seq_mat.seq,
            'indices': seq_mat.index
        }

    def generate_pre_mrna(self) -> 'Transcript':
        """Generate the pre-mRNA sequence for the transcript and store it as self.pre_mrna."""
        try:
            # Try to get sequence from FASTA if available
            seq_data = self.pull_pre_mrna_from_fasta()
            pre_mrna = SeqMat(**seq_data)
        except:
            # If no FASTA available, create empty sequence
            length = self.transcript_upper - self.transcript_lower + 1
            indices = np.arange(self.transcript_lower, self.transcript_upper + 1)
            pre_mrna = SeqMat('N' * length, indices=indices)
        
        if self.rev:
            pre_mrna.reverse_complement()
        self.pre_mrna = pre_mrna
        return self

    def generate_mature_mrna(self, inplace: bool = True) -> Union['Transcript', SeqMat]:
        """
        Generate the mature mRNA by concatenating exon regions from pre_mRNA.

        Args:
            inplace: If True, set self.mature_mrna, else return a new SeqMat

        Returns:
            The Transcript object (if inplace=True) or a SeqMat (if inplace=False)
        """
        self._fix_and_check_introns()

        if inplace:
            self.mature_mrna = self.pre_mrna.remove_regions(self.introns)
            return self

        return self.pre_mrna.remove_regions(self.introns)

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

        # Extract ORF region from mature mRNA
        if hasattr(self, 'mature_mrna') and self.mature_mrna is not None:
            # Find the positions in the mature mRNA
            orf_start_idx = np.where(self.mature_mrna.index == self.TIS)[0]
            orf_end_idx = np.where(self.mature_mrna.index == self.TTS)[0]
            
            if len(orf_start_idx) > 0 and len(orf_end_idx) > 0:
                # Extract the subsequence using the found indices
                start_idx = orf_start_idx[0]
                end_idx = orf_end_idx[0] + 1  # +1 for inclusive end
                
                # Create new SeqMat with the ORF sequence
                orf_seq = self.mature_mrna.seq[start_idx:end_idx]
                orf_indices = self.mature_mrna.index[start_idx:end_idx]
                
                return SeqMat(orf_seq, indices=orf_indices)
        
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

        if not BIOPYTHON_AVAILABLE:
            print("BioPython not available. Cannot translate to protein.")
            return self if inplace else ""

        # Translate the ORF to protein
        orf_seq = self.orf
        if isinstance(orf_seq, SeqMat):
            protein = str(Seq(orf_seq.seq).translate()).strip('*')
        else:
            protein = ""

        if inplace:
            self.protein = protein
            # Update conservation vector if available
            if hasattr(self, 'cons_vector') and self.cons_vector is not None:
                if len(self.cons_vector) != len(protein):
                    self.cons_vector = np.ones(len(protein))
            return self
        
        return protein