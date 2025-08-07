"""SeqMat - Lightning-fast genomic sequence matrix with mutation tracking"""
from __future__ import annotations

import copy
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Union, Optional, Any, ClassVar
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

from .config import get_organism_config

Mutation = Tuple[int, str, str]


def contains(array: Union[np.ndarray, List], value: Any) -> bool:
    """Check if a value is contained in an array or list"""
    if isinstance(array, np.ndarray):
        return value in array
    return value in array


@dataclass(slots=True)
class SeqMat:
    """
    Lightning-fast genomic sequence matrix with full mutation tracking,
    slicing, reverse/complement operations, and optional FASTA instantiation.

    Key features:
      - SNPs, insertions, deletions with history and sub-indexing
      - Vectorized complement & reverse-complement
      - Intuitive slicing (__getitem__) returns SeqMat clones
      - remove_regions() for excising intervals
      - Classmethod from_fasta() to load any genome FASTA
      - Per-base conservation and custom metadata support
    """
    name: str = field(default="wild_type")
    version: str = field(default="1.0")
    source: str = field(default="Unknown")
    notes: dict = field(default_factory=dict)
    
    
    # Under-the-hood storage
    seq_array: np.ndarray = field(init=False, repr=False)
    insertions: Dict[int, List[Tuple[int, str]]] = field(default_factory=lambda: defaultdict(list), init=False, repr=False)
    mutations: List[Dict[str, Any]] = field(default_factory=list, init=False, repr=False)
    mutated_positions: set[int] = field(default_factory=set, init=False, repr=False)
    rev: bool = field(default=False, init=False, repr=False)
    predicted_splicing: Optional[pd.DataFrame] = field(default=None, init=False, repr=False)

    def __init__(
        self,
        nucleotides: Union[str, np.ndarray] = "",
        indices: Optional[np.ndarray] = None,
        conservation: Optional[np.ndarray] = None,
        reference: Optional[Union[str, np.ndarray]] = None,
        name: str = 'wild_type',
        source: Optional[str] = None,
        version: str = '1.0',
        notes: Optional[dict] = None,
        rev: bool = False,
        seq: Optional[str] = None  # Alternative parameter name
    ) -> None:
        # Handle alternative parameter names
        if seq is not None and not nucleotides:
            nucleotides = seq
            
        # Metadata
        self.name = name
        self.version = version
        self.source = source or "Unknown"
        self.notes = notes or {}
        self.rev = rev
        self.predicted_splicing = None

        # Tracking
        self.insertions = defaultdict(list)
        self.mutations = []
        self.mutated_positions = set()

        # Prepare sequence
        if isinstance(nucleotides, str):
            nts = np.array(list(nucleotides), dtype='S1')
        else:
            nts = np.array(nucleotides, dtype='S1')
        L = len(nts)

        # Indices default to 1-based
        if indices is None:
            indices = np.arange(1, L+1, dtype=np.int64)
        else:
            indices = np.asarray(indices, dtype=np.int64)
        if len(indices) != L:
            raise ValueError(f"Indices length {len(indices)} != sequence length {L}")

        # Structured dtype
        dtype = np.dtype([
            ('nt', 'S1'), ('index', np.int64), ('subidx', np.int16),
            ('ref', 'S1'), ('cons', np.float32), ('mut_type', 'S10'), ('valid', bool)
        ])
        arr = np.zeros(L, dtype=dtype)
        arr['nt'] = nts
        arr['index'] = indices
        arr['subidx'] = 0

        # Reference
        if reference is None:
            arr['ref'] = nts
        else:
            if isinstance(reference, str):
                ref_arr = np.array(list(reference), dtype='S1')
            else:
                ref_arr = np.array(reference, dtype='S1')
            if len(ref_arr) != L:
                raise ValueError("Reference length mismatch")
            arr['ref'] = ref_arr

        # Conservation
        if conservation is not None:
            cons = np.asarray(conservation, dtype=np.float32)
            if len(cons) != L:
                raise ValueError("Conservation length mismatch")
            arr['cons'] = cons
        else:
            arr['cons'] = 0.0

        arr['mut_type'] = b''
        arr['valid'] = arr['nt'] != b'-'
        self.seq_array = arr
        self._refresh_mutation_state()

    @classmethod
    def from_fasta(
        cls,
        genome: str,
        chrom: str,
        start: int,
        end: int,
        **kwargs
    ) -> SeqMat:
        """
        Load a genomic interval from FASTA.
        
        Args:
            genome: Genome identifier (e.g. 'hg38')
            chrom: Chromosome name
            start: Start position (1-based)
            end: End position (1-based, inclusive)
            **kwargs: Additional arguments for SeqMat constructor
            
        Returns:
            SeqMat object containing the requested sequence
        """
        # config = get_organism_config(genome)
        
        # # Try different possible paths for FASTA files
        # chrom_source = config.get('CHROM_SOURCE')
        # fasta_path = config.get('fasta', chrom_source)
        
        # if fasta_path is None:
        #     raise ValueError(f"No FASTA path configured for genome '{genome}'")
        
        # # Check if it's a directory with individual chromosome files
        # fasta_path = Path(fasta_path)
        # if fasta_path.is_dir():
        #     # Look for individual chromosome file
        #     chrom_file = fasta_path / f"{chrom}.fasta"
        #     if not chrom_file.exists():
        #         # Try without 'chr' prefix
        #         alt_chrom = chrom.replace('chr', '')
        #         chrom_file = fasta_path / f"{alt_chrom}.fasta"
        #         if not chrom_file.exists():
        #             raise ValueError(f"Chromosome file not found: {chrom}.fasta in {fasta_path}")
        #     fasta_path = chrom_file
        # fasta = pysam.FastaFile(str(fasta_path))

        fasta = pysam.FastaFile(str('/tamir2/nicolaslynn/data/genomes/grch38/full/hg38.fa'))

        seq = fasta.fetch(f'chr{chrom}', start-1, end).upper()
        indices = np.arange(start, end+1, dtype=np.int64)
        return cls(nucleotides=seq, indices=indices, name=f"{chrom}:{start}-{end}", source=genome, **kwargs)

    @classmethod
    def from_fasta_file(cls, fasta_path: Union[str, Path], chrom: str, start: int, end: int, **kwargs) -> SeqMat:
        """Load a genomic interval directly from a FASTA file path"""
        fasta = pysam.FastaFile(str(fasta_path))
        seq = fasta.fetch(chrom, start-1, end).upper()
        indices = np.arange(start, end+1, dtype=np.int64)
        return cls(nucleotides=seq, indices=indices, name=f"{chrom}:{start}-{end}", **kwargs)

    def __len__(self) -> int:
        """Return the number of valid bases in the sequence."""
        return int(self.seq_array['valid'].sum())

    def __repr__(self) -> str:
        """Return a concise representation of the SeqMat object."""
        return f"<SeqMat {self.name}: {len(self)} bp, {len(self.mutated_positions)} muts>"

    @property
    def seq(self) -> str:
        """Return the current sequence as a string (only valid bases)."""
        return self.seq_array['nt'][self.seq_array['valid']].tobytes().decode()

    @property
    def reference_seq(self) -> str:
        """Return the reference sequence as a string (only valid bases)."""
        return self.seq_array['ref'][self.seq_array['valid']].tobytes().decode()

    @property
    def index(self) -> np.ndarray:
        """Return the genomic indices of valid bases."""
        return self.seq_array['index'][self.seq_array['valid']]

    @property
    def mutation_vector(self) -> np.ndarray:
        """Return a binary vector indicating mutated positions (1=mutated, 0=not mutated)."""
        valid_mask = self.seq_array['valid']
        indices = self.seq_array['index'][valid_mask]
        mutated_array = np.array(list(self.mutated_positions), dtype=np.int64)
        return np.isin(indices, mutated_array).astype(np.int8)

    def _refresh_mutation_state(self) -> None:
        """Update mutation tracking based on current state."""
        # Clear
        self.seq_array['mut_type'] = b''
        self.mutated_positions.clear()
        # SNPs
        snp = (self.seq_array['ref'] != self.seq_array['nt']) & self.seq_array['valid']
        self.seq_array['mut_type'][snp] = b'snp'
        self.mutated_positions.update(self.seq_array['index'][snp].tolist())
        # insertions
        for pos in self.insertions:
            self.mutated_positions.add(pos)

    def clone(self, start: Optional[int] = None, end: Optional[int] = None) -> SeqMat:
        """
        Create a copy of this SeqMat, optionally sliced to a specific range.
        
        Args:
            start: Start position (genomic coordinate)
            end: End position (genomic coordinate)
            
        Returns:
            A new SeqMat object
        """
        new = copy.copy(self)
        new.notes = copy.deepcopy(self.notes)
        new.insertions = copy.deepcopy(self.insertions)
        new.mutations = copy.deepcopy(self.mutations)
        new.mutated_positions = set(self.mutated_positions)
        if start is not None or end is not None:
            lo = start or self.index.min()
            hi = end or self.index.max()
            mask = (self.seq_array['index'] >= lo) & (self.seq_array['index'] <= hi)
            new.seq_array = self.seq_array[mask].copy()
        else:
            new.seq_array = self.seq_array.copy()
        return new

    def __getitem__(
        self, key: Union[int, slice, Tuple[int, int]]
    ) -> Union[np.void, SeqMat]:
        """
        Access sequence by genomic position or range.
        
        Args:
            key: Position (int), slice, or (start, end) tuple
            
        Returns:
            Single base record (if int) or new SeqMat (if slice/tuple)
        """
        if isinstance(key, int):
            mask = self.seq_array['index'] == key
            if not mask.any():
                raise KeyError(f"{key} not in SeqMat")
            return self.seq_array[mask][0]
        if isinstance(key, slice) or (isinstance(key, tuple) and len(key) == 2):
            if isinstance(key, slice):
                lo, hi = key.start, key.stop
            else:
                lo, hi = key
            return self.clone(lo, hi)
        raise TypeError("Index must be int, slice or (start,end)")

    def _classify_mutation(self, ref: str, alt: str) -> str:
        """Classify a mutation type based on reference and alternate alleles."""
        if ref == '-':
            return 'ins'
        if alt == '-':
            return 'del'
        if len(ref) == len(alt) == 1:
            return 'snp'
        return 'complex'

    def _validate_mutation_batch(
        self,
        muts: List[Mutation],
        *,
        allow_multiple_insertions: bool = True
    ) -> bool:
        """
        Ensure no two mutations in batch have overlapping reference spans.
        
        Args:
            muts: List of (pos, ref, alt) tuples
            allow_multiple_insertions: Whether to allow multiple insertions at same position
            
        Returns:
            True if valid, False if conflicts found
        """
        # Build a list of (start, end, idx) for each mutation
        spans = []
        for i, (pos, ref, alt) in enumerate(muts):
            if ref == '-':
                # insertion: zero-length span at pos
                start, end = pos, pos
            else:
                length = len(ref)
                start, end = pos, pos + length - 1
            spans.append((start, end, i))
        
        # Check every pair for overlap
        conflicts = []
        n = len(spans)
        for a in range(n):
            sa, ea, ia = spans[a]
            for b in range(a+1, n):
                sb, eb, ib = spans[b]
                # Overlap if intervals [sa,ea] and [sb,eb] intersect
                if not (ea < sb or eb < sa):
                    # special-case: two insertions at same pos
                    ref_a, alt_a = muts[ia][1], muts[ia][2]
                    ref_b, alt_b = muts[ib][1], muts[ib][2]
                    is_ins_a = (ref_a == '-')
                    is_ins_b = (ref_b == '-')
                    if is_ins_a and is_ins_b and allow_multiple_insertions:
                        continue
                    conflicts.append((ia, ib))
        
        if conflicts:
            lines = ["Found conflicting mutations:"]
            for ia, ib in conflicts:
                lines.append(f"  #{ia}: {muts[ia]}  <-->  #{ib}: {muts[ib]}")
            print('\n'.join(lines))
            return False
            
        return True

    def apply_mutations(
        self,
        mutations: Union[Mutation, List[Mutation]],
        *,
        permissive_ref: bool = False
    ) -> SeqMat:
        """
        Apply SNPs/insertions/deletions to the sequence.
        
        Mutations are always defined relative to the positive strand. If the sequence
        is currently on the negative strand, it will be temporarily converted to the
        positive strand for mutation application, then converted back.
        
        Args:
            mutations: Single mutation or list of (pos, ref, alt) tuples
            permissive_ref: If True, skip reference validation
            
        Returns:
            Self (for chaining)
        """
        if isinstance(mutations, tuple):
            mutations = [mutations]
        if not self._validate_mutation_batch(mutations):
            return self

        # Track original strand state
        was_on_negative_strand = self.rev
        
        # If on negative strand, convert to positive strand for mutation application
        if was_on_negative_strand:
            self.reverse_complement()  # Convert to positive strand
            
        try:
            for pos, ref, alt in mutations:
                # left normalize
                while ref and alt and ref[0] == alt[0]:
                    pos += 1
                    ref = ref[1:] or '-'
                    alt = alt[1:] or '-'
                # out-of-range
                if not contains(self.index, pos):
                    continue
                typ = self._classify_mutation(ref, alt)
                self.mutations.append({'pos': pos, 'ref': ref, 'alt': alt, 'type': typ})
                
                if typ == 'snp':
                    self._substitute(pos, ref, alt, permissive_ref)
                elif typ == 'ins':
                    self._insert(pos, alt)
                elif typ == 'del':
                    self._delete(pos, ref)
                elif typ == 'complex':
                    # first delete the reference bases
                    self._delete(pos, ref)
                    # then insert the new bases
                    self._insert(pos, alt)
        finally:
            # Always restore original strand state if we converted
            if was_on_negative_strand:
                self.reverse_complement()  # Convert back to negative strand
    
        self._refresh_mutation_state()
        return self

    def _substitute(self, pos: int, ref: str, alt: str, permissive: bool) -> None:
        """Apply a substitution mutation."""
        idx = np.where(self.seq_array['index'] == pos)[0]
        if not len(idx):
            return
        i = idx[0]
        rbytes = np.array(list(ref), dtype='S1')
        if not permissive and not np.array_equal(self.seq_array['ref'][i:i+len(rbytes)], rbytes):
            raise ValueError(f"Ref mismatch @{pos}")
        self.seq_array['nt'][i:i+len(rbytes)] = np.array(list(alt), dtype='S1')
        self.seq_array['mut_type'][i] = b'snp'

    def _insert(self, pos: int, seq: str) -> None:
        """Apply an insertion mutation."""
        subid = len(self.insertions[pos]) + 1
        self.insertions[pos].append((subid, seq))
        entries = []
        for nt in seq:
            e = np.zeros(1, dtype=self.seq_array.dtype)[0]
            e['nt'] = nt.encode()
            e['index'] = pos
            e['subidx'] = subid
            e['ref'] = b'-'
            e['cons'] = 0
            e['mut_type'] = b'ins'
            e['valid'] = True
            entries.append(e)
        i = np.searchsorted(self.seq_array['index'], pos, 'right')
        self.seq_array = np.concatenate([self.seq_array[:i], np.array(entries), self.seq_array[i:]])

    def _delete(self, pos: int, ref: str) -> None:
        """Apply a deletion mutation."""
        for i, base in enumerate(ref):
            mask = (self.seq_array['index'] == pos + i) & (self.seq_array['subidx'] == 0)
            self.seq_array['valid'][mask] = False
            self.seq_array['mut_type'][mask] = b'del'

    def complement(self, copy: bool = False) -> SeqMat:
        """
        Complement the sequence (A<->T, C<->G) in-place or return a copy.
        
        Args:
            copy: If True, return a new SeqMat object instead of modifying in-place
            
        Returns:
            Self (for chaining) or new SeqMat if copy=True
        """
        if copy:
            new = self.clone()
            # Use a temporary placeholder to avoid double replacement
            temp_array = new.seq_array['nt'].copy()
            temp_array[new.seq_array['nt'] == b'A'] = b'T'
            temp_array[new.seq_array['nt'] == b'T'] = b'A'
            temp_array[new.seq_array['nt'] == b'C'] = b'G'
            temp_array[new.seq_array['nt'] == b'G'] = b'C'
            # N and - stay the same
            new.seq_array['nt'] = temp_array
            return new
        else:
            # In-place modification
            temp_array = self.seq_array['nt'].copy()
            temp_array[self.seq_array['nt'] == b'A'] = b'T'
            temp_array[self.seq_array['nt'] == b'T'] = b'A'
            temp_array[self.seq_array['nt'] == b'C'] = b'G'
            temp_array[self.seq_array['nt'] == b'G'] = b'C'
            # N and - stay the same
            self.seq_array['nt'] = temp_array
            return self

    def reverse_complement(self, copy: bool = False) -> SeqMat:
        """
        Reverse-complement the sequence in-place or return a copy.
        
        Args:
            copy: If True, return a new SeqMat object instead of modifying in-place
            
        Returns:
            Self (for chaining) or new SeqMat if copy=True
        """
        if copy:
            new = self.clone()
            # Apply complement mapping using temporary array
            temp_array = new.seq_array['nt'].copy()
            temp_array[new.seq_array['nt'] == b'A'] = b'T'
            temp_array[new.seq_array['nt'] == b'T'] = b'A'
            temp_array[new.seq_array['nt'] == b'C'] = b'G'
            temp_array[new.seq_array['nt'] == b'G'] = b'C'
            new.seq_array['nt'] = temp_array
            
            # Reverse the sequence
            new.seq_array = new.seq_array[::-1].copy()
            new.rev = not new.rev
            return new
        else:
            # In-place modification
            # Apply complement mapping in place using temporary array
            temp_array = self.seq_array['nt'].copy()
            temp_array[self.seq_array['nt'] == b'A'] = b'T'
            temp_array[self.seq_array['nt'] == b'T'] = b'A'
            temp_array[self.seq_array['nt'] == b'C'] = b'G'
            temp_array[self.seq_array['nt'] == b'G'] = b'C'
            self.seq_array['nt'] = temp_array
            
            # Reverse the sequence
            self.seq_array = self.seq_array[::-1].copy()
            self.rev = not self.rev
            return self

    def remove_regions(self, regions: List[Tuple[int, int]]) -> SeqMat:
        """
        Excise given genomic intervals (inclusive).
        
        Args:
            regions: List of (start, end) tuples to remove
            
        Returns:
            New SeqMat with regions removed
        """
        new = self.clone()
        mask = np.ones(len(new.seq_array), bool)
        for lo, hi in regions:
            mask &= ~((new.seq_array['index'] >= min(lo, hi)) & 
                     (new.seq_array['index'] <= max(lo, hi)))
        new.seq_array = new.seq_array[mask].copy()
        return new

    def summary(self) -> str:
        """Return a summary of the SeqMat object."""
        return (f"SeqMat '{self.name}': {len(self)}bp valid, mutations={len(self.mutations)}, "
                f"inserts at {list(self.insertions.keys())}")

    def to_dict(self) -> Dict[str, Any]:
        """Convert SeqMat to a dictionary representation."""
        return {
            'name': self.name,
            'sequence': self.seq,
            'reference': self.reference_seq,
            'mutations': self.mutations,
            'length': len(self),
            'mutated_positions': list(self.mutated_positions)
        }
    
    def to_fasta(self, wrap: int = 80) -> str:
        """
        Export sequence in FASTA format.
        
        Args:
            wrap: Line length for sequence wrapping
            
        Returns:
            FASTA-formatted string
        """
        header = f">{self.name}"
        if self.mutations:
            header += f" mutations={len(self.mutations)}"
        
        seq = self.seq
        lines = [header]
        for i in range(0, len(seq), wrap):
            lines.append(seq[i:i+wrap])
        
        return '\n'.join(lines)