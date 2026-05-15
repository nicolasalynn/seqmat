"""GTF parsing and ``genes.db`` build pipeline.

Public:
- :func:`retrieve_and_parse_ensembl_annotations` — main builder used by ``seqmat setup``.
- :func:`process_transcript` — turn a transcript's GTF rows into our internal dict.
- :func:`split_fasta` — fan a multi-record FASTA out to one file per chromosome.

Reading uses pandas rather than ``gtfparse`` so we can stay on ``pyarrow>=15`` and NumPy 2.
"""
from __future__ import annotations

import logging
import pickle
import re
import sqlite3
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Tuple, Union

from tqdm import tqdm

if TYPE_CHECKING:
    import pandas as pd

_log = logging.getLogger(__name__)


def _get_pd():
    import pandas as pd
    return pd


_GTF_COLS = [
    "gene_id", "transcript_id", "feature", "start", "end", "strand", "seqname",
    "gene_biotype", "transcript_biotype", "gene_name", "tag", "protein_id",
]
_FEATURE_FILTER = {"gene", "transcript", "exon", "CDS"}
_GTF_FIXED_COLS = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
_GTF_ATTR_KEYS = ["gene_id", "transcript_id", "gene_name", "gene_biotype", "transcript_biotype", "tag", "protein_id"]


def _parse_gtf_attributes(attr_str: str) -> Dict[str, str]:
    """Parse a GTF attributes column (``key "value"; key2 "value2"; ...``) into a dict."""
    out: Dict[str, str] = {}
    if not isinstance(attr_str, str) or not attr_str.strip():
        return out
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
        out[m.group(1)] = m.group(2)
    return out


def _read_gtf_pandas(annotations_file: Union[Path, str]) -> "pd.DataFrame":
    """Read a GTF file into a DataFrame with the columns the build pipeline needs."""
    pd = _get_pd()
    path = Path(annotations_file)
    df = pd.read_csv(
        path,
        sep="\t",
        comment="#",
        header=None,
        names=_GTF_FIXED_COLS,
        dtype={"seqname": str, "source": str, "feature": str, "strand": str, "frame": str, "attributes": str},
        on_bad_lines="warn",
    )
    attr_rows = df["attributes"].map(_parse_gtf_attributes)
    for key in _GTF_ATTR_KEYS:
        df[key] = attr_rows.map(lambda d, k=key: d.get(k, ""))
    df = df.drop(columns=["attributes"])
    df["start"] = pd.to_numeric(df["start"], errors="coerce").astype("Int64")
    df["end"] = pd.to_numeric(df["end"], errors="coerce").astype("Int64")
    return df


def process_transcript(
    transcript_df: "pd.DataFrame",
    rev: bool,
    chrm: str,
    cons_data: Optional[Dict[str, Any]],
) -> Optional[Dict[str, Any]]:
    """Turn the GTF rows for one transcript into the internal transcript dict."""
    pd = _get_pd()
    if transcript_df.empty:
        return None

    transcript_rows = transcript_df[transcript_df["feature"] == "transcript"]
    if transcript_rows.empty:
        return None
    transcript = transcript_rows.iloc[0]

    exon_df = transcript_df[transcript_df["feature"] == "exon"]
    cds_df = transcript_df[transcript_df["feature"] == "CDS"]

    transcript_start, transcript_end = (
        (transcript["end"], transcript["start"]) if rev else (transcript["start"], transcript["end"])
    )

    exon_starts, exon_ends = (
        (exon_df["end"], exon_df["start"]) if rev else (exon_df["start"], exon_df["end"])
    )
    exon_starts, exon_ends = exon_starts.tolist(), exon_ends.tolist()

    if transcript_start not in exon_starts or transcript_end not in exon_ends:
        return None

    acceptors, donors = list(exon_starts), list(exon_ends)
    acceptors.remove(transcript_start)
    donors.remove(transcript_end)

    if len(acceptors) != len(donors):
        return None

    data = {
        "transcript_id": transcript["transcript_id"],
        "transcript_biotype": transcript["transcript_biotype"],
        "transcript_start": int(transcript_start),
        "transcript_end": int(transcript_end),
        "tag": transcript["tag"] if "tag" in transcript and pd.notna(transcript["tag"]) else "",
        "primary_transcript": (
            "Ensembl" in transcript["tag"] if "tag" in transcript and pd.notna(transcript["tag"]) else False
        ),
        "rev": rev,
        "chrm": chrm,
    }

    if acceptors and donors:
        data.update({"donors": donors, "acceptors": acceptors})

    if not cds_df.empty:
        cds_start, cds_end = (
            (cds_df["end"], cds_df["start"]) if rev else (cds_df["start"], cds_df["end"])
        )
        cds_start = [c for c in cds_start.tolist() if c not in acceptors]
        cds_end = [c for c in cds_end.tolist() if c not in donors]
        if len(cds_start) == 1 and len(cds_end) == 1:
            data.update({
                "TIS": cds_start[0],
                "TTS": cds_end[0],
                "protein_id": transcript["protein_id"] if "protein_id" in transcript and pd.notna(transcript["protein_id"]) else None,
            })

    if cons_data and transcript["transcript_id"] in cons_data:
        data.update({
            "cons_available": True,
            "cons_vector": cons_data[transcript["transcript_id"]]["scores"],
            "cons_seq": cons_data[transcript["transcript_id"]]["seq"],
        })
    else:
        data.update({"cons_available": False})

    return data


def _process_gene_chunk(args: Tuple[Any, Optional[Dict], Any]) -> List[Tuple[str, str, str, bytes]]:
    """Worker: process a subset of annotations (one chunk by gene_id)."""
    pd = _get_pd()
    chunk_df, cons_data, gtex_df = args
    if chunk_df.empty:
        return []
    out: List[Tuple[str, str, str, bytes]] = []
    for gene_id, gene_df in chunk_df.groupby("gene_id"):
        biotype = gene_df["gene_biotype"].unique().tolist()
        chrm = gene_df["seqname"].unique().tolist()
        strand = gene_df["strand"].unique().tolist()
        gene_attribute = gene_df[gene_df["feature"] == "gene"]
        if len(biotype) != 1 or len(chrm) != 1 or len(strand) != 1 or len(gene_attribute) != 1 or gene_attribute.empty:
            continue
        gene_attribute = gene_attribute.squeeze()
        rev = gene_attribute["strand"] == "-"
        chrm_str = str(gene_attribute["seqname"]).replace("chr", "")
        transcripts = {}
        for transcript_id, transcript_df in gene_df.groupby("transcript_id"):
            if transcript_id:
                transcript_data = process_transcript(transcript_df, rev, chrm_str, cons_data)
                if transcript_data:
                    transcripts[transcript_id] = transcript_data
        if not transcripts:
            continue
        tag = gene_attribute.get("tag")
        tag_list = tag.split(",") if pd.notna(tag) and isinstance(tag, str) else []
        gene_data = {
            "gene_name": gene_attribute["gene_name"],
            "chrm": chrm_str,
            "gene_id": gene_attribute["gene_id"],
            "gene_start": int(gene_attribute["start"]),
            "gene_end": int(gene_attribute["end"]),
            "rev": rev,
            "tag": tag_list,
            "biotype": gene_attribute["gene_biotype"],
            "transcripts": transcripts,
            "tissue_expression": {},
        }
        try:
            if gtex_df is not None and not gtex_df.empty and gene_id in gtex_df.index:
                tissue_expr = gtex_df.loc[gene_id]
                gene_data["tissue_expression"] = tissue_expr.squeeze().to_dict()
        except Exception:
            pass
        out.append((gene_data["gene_name"], gene_data["gene_id"], gene_data["biotype"], pickle.dumps(gene_data)))
    return out


def retrieve_and_parse_ensembl_annotations(
    local_path: Path,
    annotations_file: Path,
    cons_data: Optional[Dict[str, Any]],
    gtex_file: Optional[Path] = None,
    n_jobs: int = 1,
) -> None:
    """Parse Ensembl GTF annotations and create ``genes.db``.

    Set ``n_jobs > 1`` to parallelize gene-chunk processing across workers.
    The :mod:`locator` sidecar (``gene_locations.npz``) is built afterwards so
    fresh setups support position queries immediately.
    """
    pd = _get_pd()

    if gtex_file and gtex_file.exists():
        _log.info("Loading GTEx expression data...")
        gtex_df = pd.read_csv(gtex_file, delimiter="\t", header=2)
        gtex_df["Name"] = gtex_df["Name"].str.split(".").str[0]
        gtex_df = gtex_df.set_index("Name").drop(columns=["Description"])
    else:
        gtex_df = pd.DataFrame()

    _log.info("Reading GTF annotations...")
    annotations = _read_gtf_pandas(annotations_file)

    required_columns = ["gene_id", "gene_biotype", "seqname", "strand", "feature"]
    missing_columns = [col for col in required_columns if col not in annotations.columns]
    if missing_columns:
        raise ValueError(f"GTF missing columns: {missing_columns}. Available: {list(annotations.columns)}")

    annotations = annotations[annotations["feature"].isin(_FEATURE_FILTER)]
    keep_cols = [c for c in _GTF_COLS if c in annotations.columns]
    annotations = annotations[keep_cols].copy()
    _log.info("Annotations shape after filter: %s", annotations.shape)

    unique_genes = annotations["gene_id"].nunique()
    _log.info("Processing %d genes (n_jobs=%d)", unique_genes, n_jobs)

    genes_db = local_path / "genes.db"
    if genes_db.exists():
        genes_db.unlink()
    conn = sqlite3.connect(str(genes_db))
    conn.execute("CREATE TABLE genes (gene_name TEXT, gene_id TEXT, biotype TEXT, data BLOB)")
    conn.execute("CREATE INDEX idx_genes_name ON genes(gene_name)")
    conn.execute("CREATE INDEX idx_genes_id ON genes(gene_id)")
    cursor = conn.cursor()

    BATCH_SIZE = 2000
    insert_sql = "INSERT INTO genes (gene_name, gene_id, biotype, data) VALUES (?, ?, ?, ?)"

    if n_jobs > 1:
        gene_ids = annotations["gene_id"].unique()
        n_chunks = min(n_jobs, len(gene_ids))
        chunk_size = (len(gene_ids) + n_chunks - 1) // n_chunks
        chunks = [
            annotations[annotations["gene_id"].isin(gene_ids[i * chunk_size:(i + 1) * chunk_size])]
            for i in range(n_chunks)
        ]
        batch: List[Tuple[str, str, str, bytes]] = []
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            futures = {
                pool.submit(_process_gene_chunk, (chunk, cons_data, gtex_df)): i
                for i, chunk in enumerate(chunks)
            }
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Genes"):
                batch.extend(fut.result())
                if len(batch) >= BATCH_SIZE:
                    cursor.executemany(insert_sql, batch)
                    conn.commit()
                    batch.clear()
        if batch:
            cursor.executemany(insert_sql, batch)
            conn.commit()
    else:
        batch = []
        for gene_id, gene_df in tqdm(annotations.groupby("gene_id")):
            biotype = gene_df["gene_biotype"].unique().tolist()
            chrm = gene_df["seqname"].unique().tolist()
            strand = gene_df["strand"].unique().tolist()
            gene_attribute = gene_df[gene_df["feature"] == "gene"]
            if (
                len(biotype) != 1
                or len(chrm) != 1
                or len(strand) != 1
                or len(gene_attribute) != 1
                or gene_attribute.empty
            ):
                continue
            gene_attribute = gene_attribute.squeeze()
            rev = gene_attribute["strand"] == "-"
            chrm_str = str(gene_attribute["seqname"]).replace("chr", "")
            transcripts = {}
            for transcript_id, transcript_df in gene_df.groupby("transcript_id"):
                if transcript_id:
                    transcript_data = process_transcript(transcript_df, rev, chrm_str, cons_data)
                    if transcript_data:
                        transcripts[transcript_id] = transcript_data
            if not transcripts:
                continue
            tag = gene_attribute.get("tag")
            tag_list = tag.split(",") if pd.notna(tag) and isinstance(tag, str) else []
            gene_data = {
                "gene_name": gene_attribute["gene_name"],
                "chrm": chrm_str,
                "gene_id": gene_attribute["gene_id"],
                "gene_start": int(gene_attribute["start"]),
                "gene_end": int(gene_attribute["end"]),
                "rev": rev,
                "tag": tag_list,
                "biotype": gene_attribute["gene_biotype"],
                "transcripts": transcripts,
                "tissue_expression": {},
            }
            try:
                if not gtex_df.empty and gene_id in gtex_df.index:
                    gene_data["tissue_expression"] = gtex_df.loc[gene_id].squeeze().to_dict()
            except Exception:
                pass
            batch.append(
                (gene_data["gene_name"], gene_data["gene_id"], gene_data["biotype"], pickle.dumps(gene_data))
            )
            if len(batch) >= BATCH_SIZE:
                cursor.executemany(insert_sql, batch)
                conn.commit()
                batch.clear()
        if batch:
            cursor.executemany(insert_sql, batch)
            conn.commit()

    conn.close()
    _log.info("Wrote %s", genes_db)

    # Sidecar position index for chrom/position -> gene lookups.
    try:
        from .locator import _LocationIndex
        idx = _LocationIndex.build_from_db(genes_db)
        idx.save(local_path / "gene_locations.npz")
        _log.info("Wrote %s", local_path / "gene_locations.npz")
    except Exception as exc:
        _log.warning("Failed to build gene_locations.npz sidecar: %s", exc)


def split_fasta(input_file: Path, output_directory: Path, skip_existing: bool = False) -> None:
    """Split a multi-record FASTA into one ``{chromosome}.fasta`` per record.

    Skips records whose header contains an underscore (so unplaced contigs aren't written).
    """
    _log.info("Splitting FASTA file by chromosome...")
    output_directory.mkdir(exist_ok=True)

    if skip_existing:
        existing_fastas = list(output_directory.glob("*.fasta"))
        if existing_fastas:
            _log.info("Found %d existing chromosome files, skipping FASTA split", len(existing_fastas))
            return

    def _write(header: str, sequence: str) -> None:
        if not sequence or not header or "_" in header:
            return
        output_file = output_directory / f"{header}.fasta"
        if skip_existing and output_file.exists():
            _log.info("%s already exists, skipping", output_file.name)
            return
        with open(output_file, "w") as out:
            out.write(f">{header}\n{sequence}\n")
        _log.info("Created %s", output_file.name)

    with open(input_file, "r") as f:
        header = ""
        sequence_parts: List[str] = []
        for line in f:
            if line.startswith(">"):
                _write(header, "".join(sequence_parts))
                header = line[1:].split()[0]
                sequence_parts = []
            else:
                sequence_parts.append(line.strip())
        _write(header, "".join(sequence_parts))


__all__ = [
    "retrieve_and_parse_ensembl_annotations",
    "process_transcript",
    "split_fasta",
]
