"""Gene class for representing genomic genes with associated transcripts"""
import logging
from collections import Counter
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple, Union

from ._io import unload_pickle
from .config import get_default_organism, get_organism_config
from .errors import GeneNotFoundError, OrganismNotConfiguredError
from .locator import PosArg, gene_names_at_position
from .transcript import Transcript

_log = logging.getLogger(__name__)


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

    def __iter__(self) -> Iterator[Transcript]:
        """Allow iteration over the gene's transcripts, yielding Transcript objects."""
        for tid, annotations in self.transcripts.items():
            yield Transcript(annotations, organism=self.organism)

    def __getitem__(self, item: str) -> Optional[Transcript]:
        """Get a transcript by ID. Returns ``None`` if the ID isn't an annotated transcript."""
        if item not in self.transcripts:
            _log.warning("%s is not an annotated transcript of %s", item, self.gene_name)
            return None
        return Transcript(self.transcripts[item], organism=self.organism)

    @classmethod
    def _load_dict(cls, gene_name: str, organism: str) -> Optional[Dict[str, Any]]:
        """Load the raw gene dict from LMDB → SQLite → per-gene pickle, in that order.

        Returns ``None`` when the gene isn't found in any backend.
        """
        try:
            from .lmdb_store import load_gene_from_lmdb
            data = load_gene_from_lmdb(gene_name, organism)
            if data is not None:
                return data
        except ImportError:
            pass
        try:
            from .sqlite_store import load_gene_from_sqlite
            data = load_gene_from_sqlite(gene_name, organism)
            if data is not None:
                return data
        except ImportError:
            pass
        try:
            config = get_organism_config(organism)
        except ValueError as exc:
            raise OrganismNotConfiguredError(
                f"Organism '{organism}' not configured. Run setup_genomics_data() first."
            ) from exc
        mrna_path = Path(config["MRNA_PATH"])
        if mrna_path.exists():
            for biotype_dir in mrna_path.iterdir():
                if biotype_dir.is_dir():
                    for pkl in biotype_dir.glob(f"*_{gene_name}.pkl"):
                        return unload_pickle(pkl)
        return None

    @classmethod
    def get(cls, gene_name: str, organism: Optional[str] = None) -> 'Gene':
        """Load a gene, raising on failure.

        Raises:
            OrganismNotConfiguredError: when the organism has no data set up.
            GeneNotFoundError: when the organism is configured but the gene is missing.
        """
        if organism is None:
            organism = get_default_organism()
        data = cls._load_dict(gene_name, organism)
        if data is None:
            raise GeneNotFoundError(
                f"Gene '{gene_name}' not found in organism '{organism}'."
            )
        return cls(
            gene_name=data.get('gene_name'),
            gene_id=data.get('gene_id'),
            rev=data.get('rev'),
            chrm=data.get('chrm'),
            transcripts=data.get('transcripts', {}),
            organism=organism,
        )

    @classmethod
    def from_file(cls, gene_name: str, organism: Optional[str] = None) -> Optional['Gene']:
        """Load a gene by name, returning ``None`` if not found.

        Load order: LMDB (if installed) → SQLite (``genes.db``) → per-gene pickle files.
        Prefer :meth:`get` if you'd rather have a typed exception than ``None``.
        """
        if organism is None:
            organism = get_default_organism()
        try:
            return cls.get(gene_name, organism=organism)
        except OrganismNotConfiguredError as exc:
            _log.warning("%s", exc)
            return None
        except GeneNotFoundError as exc:
            _log.warning("%s", exc)
            return None

    @classmethod
    def from_position(
        cls,
        chrm: str,
        pos: PosArg,
        organism: Optional[str] = None,
    ) -> List['Gene']:
        """Return all genes overlapping a point or range on a chromosome.

        Args:
            chrm: Chromosome (e.g. "12", "chr12", "X"). Leading 'chr' is stripped.
            pos: Either an int position or a (start, end) tuple (inclusive).
            organism: Organism build (uses default if None).

        Returns:
            List of Gene objects, possibly empty. Returned in ascending start order.
        """
        if organism is None:
            organism = get_default_organism()
        names = gene_names_at_position(chrm, pos, organism=organism)
        genes: List['Gene'] = []
        for name in names:
            g = cls.from_file(name, organism=organism)
            if g is not None:
                genes.append(g)
        return genes

    def splice_sites(self) -> Tuple[Counter, Counter]:
        """Return (Counter of acceptors, Counter of donors) across all transcripts."""
        acceptors: List[Any] = []
        donors: List[Any] = []
        for transcript in self.transcripts.values():
            acceptors.extend(transcript.get('acceptors', []))
            donors.extend(transcript.get('donors', []))

        return Counter(acceptors), Counter(donors)

    def transcript(self, tid: Optional[str] = None) -> Optional[Transcript]:
        """Return Transcript by ID, or primary transcript if tid is None."""
        if tid is None:
            tid = self.primary_transcript
            
        if tid is None or tid not in self.transcripts:
            return None

        return Transcript(self.transcripts[tid], organism=self.organism)

    @property
    def primary_transcript(self) -> Optional[str]:
        """Primary transcript ID, or first protein-coding transcript, or None."""
        if hasattr(self, "_primary_transcript"):
            return self._primary_transcript
        primary = [k for k, v in self.transcripts.items() if v.get("primary_transcript")]
        if primary:
            self._primary_transcript = primary[0]
            return self._primary_transcript
        protein_coding = [k for k, v in self.transcripts.items()
                         if v.get("transcript_biotype") == "protein_coding"]
        if protein_coding:
            self._primary_transcript = protein_coding[0]
            return self._primary_transcript
        self._primary_transcript = None
        return None