"""Gene class for representing genomic genes with associated transcripts"""
import copy
from typing import Any, Dict, List, Tuple, Optional, Iterator, Union
from collections import Counter
from pathlib import Path

from .config import get_organism_config, get_default_organism
from .utils import unload_pickle
from .transcript import Transcript


class Gene:
    """
    A class representing a Gene, with associated transcripts and metadata.

    Attributes:
        organism (str): The organism build (e.g. 'hg38').
        transcripts (dict): A dictionary of transcript annotations keyed by transcript ID.
        gene_name (str): The name of the gene.
        gene_id (str): The unique identifier for the gene.
        chrm (str): The chromosome on which the gene resides.
        rev (bool): Whether the gene is on the reverse strand.
    """

    def __init__(self, gene_name: str, gene_id: str, rev: bool, chrm: str, 
                 transcripts: Optional[Dict[str, Any]] = None, organism: Optional[str] = None):
        """
        Initialize a Gene instance.

        Args:
            gene_name: Name of the gene
            gene_id: Unique identifier for the gene
            rev: Whether gene is on reverse strand
            chrm: Chromosome identifier
            transcripts: Dictionary of transcript annotations
            organism: Organism reference build (default from config)
        """
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.rev = rev
        self.chrm = chrm
        self.organism = organism if organism is not None else get_default_organism()
        self.transcripts = transcripts if transcripts is not None else {}

    def __repr__(self) -> str:
        """Official string representation of the Gene object."""
        return f"Gene({self.gene_name})"

    def __str__(self) -> str:
        """User-friendly string representation of the Gene object."""
        return f"Gene: {self.gene_name}, ID: {self.gene_id}, Chr: {self.chrm}, Transcripts: {len(self.transcripts)}"

    def __len__(self) -> int:
        """Returns the number of transcripts associated with this gene."""
        return len(self.transcripts)

    def __copy__(self):
        """Returns a shallow copy of the Gene object."""
        return copy.copy(self)

    def __deepcopy__(self, memo):
        """Returns a deep copy of the Gene object."""
        return copy.deepcopy(self, memo)

    def __iter__(self) -> Iterator[Transcript]:
        """Allow iteration over the gene's transcripts, yielding Transcript objects."""
        for tid, annotations in self.transcripts.items():
            yield Transcript(annotations, organism=self.organism)

    def __getitem__(self, item: str) -> Optional[Transcript]:
        """Get a transcript by ID."""
        if item not in self.transcripts:
            print(f"{item} not an annotated transcript of this gene.")
            return None
        return Transcript(self.transcripts[item], organism=self.organism)

    @classmethod
    def from_file(cls, gene_name: str, organism: Optional[str] = None) -> Optional['Gene']:
        """
        Load gene data from file.
        
        Args:
            gene_name: Name of the gene to load
            organism: Organism reference build
            
        Returns:
            Gene object or None if not found
        """
        if organism is None:
            organism = get_default_organism()
        try:
            config = get_organism_config(organism)
        except ValueError:
            print(f"Organism '{organism}' not configured. Run setup_genomics_data() first.")
            return None
            
        # Search through all biotype folders in the configured organism MRNA path
        mrna_path = Path(config['MRNA_PATH'])
        gene_files = []
        
        # Look through all biotype subdirectories
        if mrna_path.exists():
            for biotype_dir in mrna_path.iterdir():
                if biotype_dir.is_dir():
                    # Search for gene files matching the name
                    matching_files = list(biotype_dir.glob(f'*_{gene_name}.pkl'))
                    gene_files.extend(matching_files)
        
        if not gene_files:
            print(f"No files available for gene '{gene_name}'.")
            return None

        # Load gene data from the first matching file
        data = unload_pickle(gene_files[0])
        
        return cls(
            gene_name=data.get('gene_name'),
            gene_id=data.get('gene_id'),
            rev=data.get('rev'),
            chrm=data.get('chrm'),
            transcripts=data.get('transcripts', {}),
            organism=organism
        )

    def splice_sites(self) -> Tuple[Counter, Counter]:
        """
        Aggregates splice sites (acceptors and donors) from all transcripts.

        Returns:
            tuple(Counter, Counter): A tuple of two Counters for acceptors and donors.
        """
        acceptors: List[Any] = []
        donors: List[Any] = []

        # Collect acceptor and donor sites from each transcript
        for transcript in self.transcripts.values():
            acceptors.extend(transcript.get('acceptors', []))
            donors.extend(transcript.get('donors', []))

        return Counter(acceptors), Counter(donors)

    def transcript(self, tid: Optional[str] = None) -> Optional[Transcript]:
        """
        Retrieve a Transcript object by ID, or the primary transcript if no ID is given.

        Args:
            tid: Transcript ID. If None, returns primary transcript.

        Returns:
            The Transcript object with the given ID or the primary transcript.
        """
        if tid is None:
            tid = self.primary_transcript
            
        if tid is None or tid not in self.transcripts:
            return None

        return Transcript(self.transcripts[tid], organism=self.organism)

    @property
    def primary_transcript(self) -> Optional[str]:
        """
        Returns the primary transcript ID for this gene.
        
        Returns:
            The primary transcript ID or None if not available.
        """
        # If already calculated, return it
        if hasattr(self, '_primary_transcript'):
            return self._primary_transcript

        # Try to find a primary transcript
        primary_transcripts = [k for k, v in self.transcripts.items() 
                             if v.get('primary_transcript')]
        if primary_transcripts:
            self._primary_transcript = primary_transcripts[0]
            return self._primary_transcript

        # Fallback: find a protein-coding transcript
        protein_coding = [k for k, v in self.transcripts.items() 
                        if v.get('transcript_biotype') == 'protein_coding']
        if protein_coding:
            self._primary_transcript = protein_coding[0]
            return self._primary_transcript

        # No primary or protein-coding transcript found
        self._primary_transcript = None
        return None