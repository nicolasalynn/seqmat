"""Utility functions for SeqMat package"""
from __future__ import annotations

import os
import pickle
import json
import re
import sqlite3
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
import numpy as np
from pathlib import Path
from typing import TYPE_CHECKING, Any, Union, List, Optional, Dict, Tuple
from tqdm import tqdm

if TYPE_CHECKING:
    import pandas as pd
else:
    def _get_pd():
        import pandas as pd
        return pd
import requests
import gzip
import shutil

# GTF reading: we use a pandas-based reader so we don't depend on gtfparse (which pins pyarrow<14.1).
# That allows seqmat to require pyarrow>=15 and NumPy 2.

try:
    import polars as pl
    POLARS_AVAILABLE = True
except ImportError:
    POLARS_AVAILABLE = False

from .config import (
    save_config,
    load_config,
    get_default_organism,
    get_directory_config,
    get_available_organisms,
    DEFAULT_ORGANISM_DATA,
    get_organism_info,
    get_data_dir,
    get_prebuilt_data_base_url,
)

# Columns needed for GTF processing (drop the rest to save memory). Missing cols handled in code.
_GTF_COLS = [
    'gene_id', 'transcript_id', 'feature', 'start', 'end', 'strand', 'seqname',
    'gene_biotype', 'transcript_biotype', 'gene_name', 'tag', 'protein_id'
]
_FEATURE_FILTER = {'gene', 'transcript', 'exon', 'CDS'}

_GTF_FIXED_COLS = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
_GTF_ATTR_KEYS = ['gene_id', 'transcript_id', 'gene_name', 'gene_biotype', 'transcript_biotype', 'tag', 'protein_id']


def _parse_gtf_attributes(attr_str: str) -> Dict[str, str]:
    """Parse GTF attributes column (key "value"; key2 "value2"; ...) into a dict."""
    out: Dict[str, str] = {}
    if not isinstance(attr_str, str) or not attr_str.strip():
        return out
    for m in re.finditer(r'(\w+)\s+"([^"]*)"', attr_str):
        out[m.group(1)] = m.group(2)
    return out


def _read_gtf_pandas(annotations_file: Union[Path, str]) -> "pd.DataFrame":
    """Read a GTF file into a pandas DataFrame with columns needed for annotation processing.
    Does not require gtfparse, so we can use pyarrow>=15 and NumPy 2.
    """
    pd = _get_pd()
    path = Path(annotations_file)
    df = pd.read_csv(
        path,
        sep='\t',
        comment='#',
        header=None,
        names=_GTF_FIXED_COLS,
        dtype={'seqname': str, 'source': str, 'feature': str, 'strand': str, 'frame': str, 'attributes': str},
        on_bad_lines='warn',
    )
    # Parse attributes into columns
    attr_rows = df['attributes'].map(_parse_gtf_attributes)
    for key in _GTF_ATTR_KEYS:
        df[key] = attr_rows.map(lambda d, k=key: d.get(k, ''))
    df = df.drop(columns=['attributes'])
    df['start'] = pd.to_numeric(df['start'], errors='coerce').astype('Int64')
    df['end'] = pd.to_numeric(df['end'], errors='coerce').astype('Int64')
    return df


def _process_gene_chunk(args: Tuple[Any, Optional[Dict], Any]) -> List[Tuple[str, str, str, bytes]]:
    """Process a subset of annotations (one chunk by gene_id). Returns list of (gene_name, gene_id, biotype, blob)."""
    pd = _get_pd()
    chunk_df, cons_data, gtex_df = args
    if chunk_df.empty:
        return []
    out: List[Tuple[str, str, str, bytes]] = []
    for gene_id, gene_df in chunk_df.groupby('gene_id'):
        biotype = gene_df['gene_biotype'].unique().tolist()
        chrm = gene_df['seqname'].unique().tolist()
        strand = gene_df['strand'].unique().tolist()
        gene_attribute = gene_df[gene_df['feature'] == 'gene']
        if len(biotype) != 1 or len(chrm) != 1 or len(strand) != 1 or len(gene_attribute) != 1 or gene_attribute.empty:
            continue
        gene_attribute = gene_attribute.squeeze()
        rev = gene_attribute['strand'] == '-'
        chrm_str = str(gene_attribute['seqname']).replace('chr', '')
        transcripts = {}
        for transcript_id, transcript_df in gene_df.groupby('transcript_id'):
            if transcript_id:
                transcript_data = process_transcript(transcript_df, rev, chrm_str, cons_data)
                if transcript_data:
                    transcripts[transcript_id] = transcript_data
        if not transcripts:
            continue
        tag = gene_attribute.get('tag')
        tag_list = tag.split(',') if pd.notna(tag) and isinstance(tag, str) else []
        gene_data = {
            'gene_name': gene_attribute['gene_name'],
            'chrm': chrm_str,
            'gene_id': gene_attribute['gene_id'],
            'gene_start': int(gene_attribute['start']),
            'gene_end': int(gene_attribute['end']),
            'rev': rev,
            'tag': tag_list,
            'biotype': gene_attribute['gene_biotype'],
            'transcripts': transcripts,
            'tissue_expression': {},
        }
        try:
            if gtex_df is not None and not gtex_df.empty and gene_id in gtex_df.index:
                tissue_expr = gtex_df.loc[gene_id]
                gene_data['tissue_expression'] = tissue_expr.squeeze().to_dict()
        except Exception:
            pass
        out.append((gene_data['gene_name'], gene_data['gene_id'], gene_data['biotype'], pickle.dumps(gene_data)))
    return out


def dump_pickle(path: Union[str, Path], data: Any) -> None:
    """Save data to a pickle file"""
    with open(path, 'wb') as f:
        pickle.dump(data, f)


def unload_pickle(path: Union[str, Path]) -> Any:
    """Load data from a pickle file"""
    with open(path, 'rb') as f:
        return pickle.load(f)


def load_conservation(path: Union[str, Path]) -> Dict[str, Any]:
    """
    Load conservation data from .db (SQLite) or .pkl.
    Returns dict[transcript_id] -> {'scores': np.ndarray, 'seq': str/bytes}.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".db":
        conn = sqlite3.connect(str(path))
        conn.row_factory = sqlite3.Row
        cursor = conn.execute("SELECT transcript_id, data FROM conservation")
        out: Dict[str, Any] = {}
        for row in cursor:
            out[str(row["transcript_id"])] = pickle.loads(row["data"])
        conn.close()
        return out
    data = unload_pickle(path)
    if not isinstance(data, dict):
        raise ValueError("Conservation pickle must be a dict[transcript_id] -> {scores, seq}")
    return data


def convert_conservation_pkl_to_db(pkl_path: Union[str, Path], db_path: Optional[Union[str, Path]] = None) -> Path:
    """
    Convert conservation.pkl to conservation.db (SQLite, same pattern as genes.db).
    Returns path to the written .db file.
    """
    pkl_path = Path(pkl_path)
    db_path = Path(db_path) if db_path is not None else pkl_path.with_suffix(".db")
    if db_path.exists():
        db_path.unlink()
    data = unload_pickle(pkl_path)
    if not isinstance(data, dict):
        raise ValueError("Conservation pickle must be a dict[transcript_id] -> {scores, seq}")
    conn = sqlite3.connect(str(db_path))
    conn.execute("CREATE TABLE conservation (transcript_id TEXT PRIMARY KEY, data BLOB)")
    for transcript_id, blob in data.items():
        conn.execute("INSERT INTO conservation (transcript_id, data) VALUES (?, ?)", (str(transcript_id), pickle.dumps(blob)))
    conn.commit()
    conn.close()
    return db_path


def dump_json(path: Union[str, Path], data: Any) -> None:
    """Save data to a JSON file"""
    with open(path, 'w') as f:
        json.dump(data, f, indent=2)


def unload_json(path: Union[str, Path]) -> Any:
    """Load data from a JSON file"""
    with open(path, 'r') as f:
        return json.load(f)


def download(external_url: str, local_path: Path, skip_existing: bool = False) -> Path:
    """Download a file from URL to local path."""
    local_file = Path(external_url).name
    local_file_path = Path(local_path) / local_file
    
    # Check if file already exists
    if skip_existing and local_file_path.exists():
        print(f"✓ File already exists: {local_file_path.name}")
        return local_file_path
    
    print(f"Downloading {external_url}")
    
    try:
        response = requests.get(external_url, stream=True)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        block_size = 8192
        
        with open(local_file_path, 'wb') as f:
            with tqdm(total=total_size, unit='iB', unit_scale=True) as pbar:
                for chunk in response.iter_content(block_size):
                    if chunk:
                        f.write(chunk)
                        pbar.update(len(chunk))
    except Exception as e:
        print(f"Error during download: {e}")
        raise

    return local_file_path


def download_and_ungzip(external_url: str, local_path: Path, skip_existing: bool = False) -> Path:
    """Download and decompress a gzipped file."""
    # Check if uncompressed file already exists
    local_file = Path(external_url).name
    output_file = Path(local_path) / local_file.rstrip('.gz')
    
    if skip_existing and output_file.exists():
        print(f"✓ Uncompressed file already exists: {output_file.name}")
        return output_file
    
    gzipped_file = download(external_url, local_path, skip_existing=skip_existing)
    
    # Check if we need to decompress
    if gzipped_file.suffix == '.gz' and not output_file.exists():
        print(f"Decompressing {gzipped_file.name}")
        
        with gzip.open(gzipped_file, 'rb') as f_in:
            with open(output_file, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        
        # Remove the gzipped file
        gzipped_file.unlink()
    
    return output_file


def process_transcript(transcript_df: "pd.DataFrame", rev: bool, chrm: str,
                      cons_data: Optional[Dict[str, Any]]) -> Optional[Dict[str, Any]]:
    """Process transcript data from GTF dataframe."""
    pd = _get_pd()
    if transcript_df.empty:
        return None

    transcript_rows = transcript_df[transcript_df['feature'] == 'transcript']
    if transcript_rows.empty:
        return None
    transcript = transcript_rows.iloc[0]

    exon_df = transcript_df[transcript_df['feature'] == 'exon']
    cds_df = transcript_df[transcript_df['feature'] == 'CDS']

    # Simplifying start and end assignments
    transcript_start, transcript_end = (transcript['end'], transcript['start']) if rev else (transcript['start'], transcript['end'])

    # Handling exons
    exon_starts, exon_ends = (exon_df['end'], exon_df['start']) if rev else (exon_df['start'], exon_df['end'])
    exon_starts, exon_ends = exon_starts.tolist(), exon_ends.tolist()

    if transcript_start not in exon_starts or transcript_end not in exon_ends:
        return None

    acceptors, donors = list(exon_starts), list(exon_ends)
    acceptors.remove(transcript_start)
    donors.remove(transcript_end)

    if len(acceptors) != len(donors):
        return None

    data = {
        'transcript_id': transcript['transcript_id'],
        'transcript_biotype': transcript['transcript_biotype'],
        'transcript_start': int(transcript_start),
        'transcript_end': int(transcript_end),
        'tag': transcript['tag'] if 'tag' in transcript and pd.notna(transcript['tag']) else '',
        'primary_transcript': 'Ensembl' in transcript['tag'] if 'tag' in transcript and pd.notna(transcript['tag']) else False,
        'rev': rev,
        'chrm': chrm
    }

    if acceptors and donors:
        data.update({'donors': donors, 'acceptors': acceptors})

    # Handling CDS
    if not cds_df.empty:
        cds_start, cds_end = (cds_df['end'], cds_df['start']) if rev else (cds_df['start'], cds_df['end'])
        cds_start = [c for c in cds_start.tolist() if c not in acceptors]
        cds_end = [c for c in cds_end.tolist() if c not in donors]

        if len(cds_start) == 1 and len(cds_end) == 1:
            data.update({
                'TIS': cds_start[0],
                'TTS': cds_end[0],
                'protein_id': transcript['protein_id'] if 'protein_id' in transcript and pd.notna(transcript['protein_id']) else None
            })

    # Add conservation data if available
    if cons_data and transcript['transcript_id'] in cons_data:
        data.update({
            'cons_available': True,
            'cons_vector': cons_data[transcript['transcript_id']]['scores'],
            'cons_seq': cons_data[transcript['transcript_id']]['seq']
        })
    else:
        data.update({'cons_available': False})

    return data


def retrieve_and_parse_ensembl_annotations(
    local_path: Path,
    annotations_file: Path,
    cons_data: Optional[Dict[str, Any]],
    gtex_file: Optional[Path] = None,
    n_jobs: int = 1,
) -> None:
    """Parse Ensembl GTF annotations and create genes.db. Use n_jobs > 1 to parallelize (faster)."""
    pd = _get_pd()

    # Load GTEx expression data if available
    if gtex_file and gtex_file.exists():
        print("Loading GTEx expression data...")
        gtex_df = pd.read_csv(gtex_file, delimiter='\t', header=2)
        gtex_df['Name'] = gtex_df['Name'].str.split('.').str[0]
        gtex_df = gtex_df.set_index('Name').drop(columns=['Description'])
    else:
        gtex_df = pd.DataFrame()

    print("Reading GTF annotations...")
    annotations = _read_gtf_pandas(annotations_file)

    required_columns = ['gene_id', 'gene_biotype', 'seqname', 'strand', 'feature']
    missing_columns = [col for col in required_columns if col not in annotations.columns]
    if missing_columns:
        raise ValueError(f"GTF missing columns: {missing_columns}. Available: {list(annotations.columns)}")

    # Pre-filter: keep only gene/transcript/exon/CDS and needed columns (speeds up groupby and iteration)
    annotations = annotations[annotations['feature'].isin(_FEATURE_FILTER)]
    keep_cols = [c for c in _GTF_COLS if c in annotations.columns]
    annotations = annotations[keep_cols].copy()
    print(f"Annotations shape after filter: {annotations.shape}")

    unique_genes = annotations['gene_id'].nunique()
    print(f"Processing {unique_genes} genes (n_jobs={n_jobs})...")

    genes_db = local_path / 'genes.db'
    if genes_db.exists():
        genes_db.unlink()
    conn = sqlite3.connect(str(genes_db))
    conn.execute("CREATE TABLE genes (gene_name TEXT, gene_id TEXT, biotype TEXT, data BLOB)")
    conn.execute("CREATE INDEX idx_genes_name ON genes(gene_name)")
    conn.execute("CREATE INDEX idx_genes_id ON genes(gene_id)")
    cursor = conn.cursor()

    BATCH_SIZE = 2000

    if n_jobs > 1:
        # Parallel: split by gene_id, process chunks in worker processes
        gene_ids = annotations['gene_id'].unique()
        n_chunks = min(n_jobs, len(gene_ids))
        chunk_size = (len(gene_ids) + n_chunks - 1) // n_chunks
        chunks = [
            annotations[annotations['gene_id'].isin(gene_ids[i * chunk_size:(i + 1) * chunk_size])]
            for i in range(n_chunks)
        ]
        batch: List[Tuple[str, str, str, bytes]] = []
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            futures = {
                pool.submit(_process_gene_chunk, (chunk, cons_data, gtex_df)): i
                for i, chunk in enumerate(chunks)
            }
            for fut in tqdm(as_completed(futures), total=len(futures), desc="Genes"):
                rows = fut.result()
                batch.extend(rows)
                if len(batch) >= BATCH_SIZE:
                    cursor.executemany(
                        "INSERT INTO genes (gene_name, gene_id, biotype, data) VALUES (?, ?, ?, ?)",
                        batch,
                    )
                    conn.commit()
                    batch.clear()
        if batch:
            cursor.executemany(
                "INSERT INTO genes (gene_name, gene_id, biotype, data) VALUES (?, ?, ?, ?)",
                batch,
            )
            conn.commit()
    else:
        # Serial path with batched INSERT
        batch: List[Tuple[str, str, str, bytes]] = []
        for gene_id, gene_df in tqdm(annotations.groupby('gene_id')):
            biotype = gene_df['gene_biotype'].unique().tolist()
            chrm = gene_df['seqname'].unique().tolist()
            strand = gene_df['strand'].unique().tolist()
            gene_attribute = gene_df[gene_df['feature'] == 'gene']
            if len(biotype) != 1 or len(chrm) != 1 or len(strand) != 1 or len(gene_attribute) != 1 or gene_attribute.empty:
                continue
            gene_attribute = gene_attribute.squeeze()
            rev = gene_attribute['strand'] == '-'
            chrm_str = str(gene_attribute['seqname']).replace('chr', '')
            transcripts = {}
            for transcript_id, transcript_df in gene_df.groupby('transcript_id'):
                if transcript_id:
                    transcript_data = process_transcript(transcript_df, rev, chrm_str, cons_data)
                    if transcript_data:
                        transcripts[transcript_id] = transcript_data
            if not transcripts:
                continue
            tag = gene_attribute.get('tag')
            tag_list = tag.split(',') if pd.notna(tag) and isinstance(tag, str) else []
            gene_data = {
                'gene_name': gene_attribute['gene_name'],
                'chrm': chrm_str,
                'gene_id': gene_attribute['gene_id'],
                'gene_start': int(gene_attribute['start']),
                'gene_end': int(gene_attribute['end']),
                'rev': rev,
                'tag': tag_list,
                'biotype': gene_attribute['gene_biotype'],
                'transcripts': transcripts,
                'tissue_expression': {},
            }
            try:
                if not gtex_df.empty and gene_id in gtex_df.index:
                    gene_data['tissue_expression'] = gtex_df.loc[gene_id].squeeze().to_dict()
            except Exception:
                pass
            batch.append((gene_data['gene_name'], gene_data['gene_id'], gene_data['biotype'], pickle.dumps(gene_data)))
            if len(batch) >= BATCH_SIZE:
                cursor.executemany(
                    "INSERT INTO genes (gene_name, gene_id, biotype, data) VALUES (?, ?, ?, ?)",
                    batch,
                )
                conn.commit()
                batch.clear()
        if batch:
            cursor.executemany(
                "INSERT INTO genes (gene_name, gene_id, biotype, data) VALUES (?, ?, ?, ?)",
                batch,
            )
            conn.commit()

    conn.close()
    print(f"Wrote {genes_db}")

def split_fasta(input_file: Path, output_directory: Path, skip_existing: bool = False) -> None:
    """Split a FASTA file into individual chromosome files."""
    print("Splitting FASTA file by chromosome...")
    
    output_directory.mkdir(exist_ok=True)
    
    # Check if chromosomes already exist
    if skip_existing:
        existing_fastas = list(output_directory.glob("*.fasta"))
        if existing_fastas:
            print(f"✓ Found {len(existing_fastas)} existing chromosome files, skipping FASTA split")
            return
    
    with open(input_file, 'r') as file:
        sequence = ''
        header = ''
        
        for line in file:
            if line.startswith('>'):
                if sequence and header and '_' not in header:
                    # Write previous sequence
                    output_file = output_directory / f"{header}.fasta"
                    if skip_existing and output_file.exists():
                        print(f"✓ {output_file.name} already exists, skipping")
                    else:
                        with open(output_file, 'w') as out:
                            out.write(f'>{header}\n{sequence}\n')
                        print(f"Created {output_file.name}")
                    
                # Start new sequence
                header = line[1:].split()[0]
                sequence = ''
            else:
                sequence += line.strip()
        
        # Write last sequence
        if sequence and header and '_' not in header:
            output_file = output_directory / f"{header}.fasta"
            with open(output_file, 'w') as out:
                out.write(f'>{header}\n{sequence}\n')
            print(f"Created {output_file.name}")


def download_genome_data(organism: str, base_path: Path, skip_existing: bool = False) -> Dict[str, Path]:
    """Download all required data files for a specific organism."""
    # Always get URLs from DEFAULT_ORGANISM_DATA, not from config
    from .config import DEFAULT_ORGANISM_DATA
    
    if organism not in DEFAULT_ORGANISM_DATA:
        raise ValueError(f"Organism {organism} not supported for download. Available: {list(DEFAULT_ORGANISM_DATA.keys())}")
    
    organism_data = DEFAULT_ORGANISM_DATA[organism]
    urls = organism_data.get('urls', {})
    
    if not urls:
        raise ValueError(f"No URLs configured for organism {organism}")
    
    # Map config URLs to expected keys (conservation comes from prebuilt bucket below)
    url_mapping = {
        'fasta': 'fasta_url',
        'gtf': 'ensembl_url',
        'gtex': 'expression_url'
    }
    
    # Convert config format to expected format
    mapped_urls = {}
    for config_key, url in urls.items():
        if config_key in url_mapping:
            mapped_urls[url_mapping[config_key]] = url
    
    if not mapped_urls:
        raise ValueError(f"No valid URLs found for organism {organism}. Available keys: {list(urls.keys())}")
    
    urls = mapped_urls
    files = {}
    
    # Download Ensembl annotations
    if not urls.get('ensembl_url'):
        raise ValueError(f"No GTF URL configured for organism {organism}")
    files['ensembl_file'] = download_and_ungzip(urls['ensembl_url'], base_path, skip_existing=skip_existing)

    # Conservation: from prebuilt bucket (same as default download), optional
    files['cons_file'] = None
    base_url = get_prebuilt_data_base_url()
    for ext in (".db", ".pkl"):
        try:
            files['cons_file'] = download(
                f"{base_url}/{organism}/conservation{ext}", base_path, skip_existing=skip_existing
            )
            break
        except Exception:
            continue
    
    # Download expression data if available
    if urls.get('expression_url'):
        files['gtex_file'] = download_and_ungzip(urls['expression_url'], base_path, skip_existing=skip_existing)
    else:
        files['gtex_file'] = None
    
    # Download genome FASTA
    if not urls.get('fasta_url'):
        raise ValueError(f"No FASTA URL configured for organism {organism}")
    files['fasta_file'] = download_and_ungzip(urls['fasta_url'], base_path, skip_existing=skip_existing)
    
    return files


class PrebuiltDataUnavailableError(Exception):
    """Raised when prebuilt data cannot be downloaded (e.g. 404). Use build-from-sources instead."""


def download_prebuilt_data(
    organism: str,
    base_path: Path,
    skip_existing: bool = False,
) -> Dict[str, Path]:
    """
    Download prebuilt genes.db and FASTA from the SeqMat S3 bucket.
    Layout: {base}/{organism}/genes.db, {base}/{organism}/{organism}.fa.gz.
    Optionally {base}/{organism}/conservation.db or conservation.pkl if present.
    """
    from .config import DEFAULT_ORGANISM_DATA

    if organism not in DEFAULT_ORGANISM_DATA:
        raise ValueError(
            f"Organism {organism} not supported. Available: {list(DEFAULT_ORGANISM_DATA.keys())}"
        )
    base_path = Path(base_path)
    base_path.mkdir(parents=True, exist_ok=True)
    base_url = get_prebuilt_data_base_url()
    genes_db_url = f"{base_url}/{organism}/genes.db"
    fasta_gz_url = f"{base_url}/{organism}/{organism}.fa.gz"
    files: Dict[str, Path] = {}

    def _download_or_raise(url: str, dest_dir: Path, is_gz: bool = False) -> Path:
        try:
            if is_gz:
                return download_and_ungzip(url, dest_dir, skip_existing=skip_existing)
            return download(url, dest_dir, skip_existing=skip_existing)
        except requests.HTTPError as e:
            if e.response is not None and e.response.status_code == 404:
                raise PrebuiltDataUnavailableError(
                    f"Prebuilt data not found at {url}. "
                    "Upload genes.db and FASTA to the S3 bucket, or run setup with "
                    "--build-from-sources to generate data from GTF/FASTA sources."
                ) from e
            raise
        except Exception as e:
            if "404" in str(e).lower() or "not found" in str(e).lower():
                raise PrebuiltDataUnavailableError(
                    f"Prebuilt data not found for {organism}. "
                    "Run setup with --build-from-sources to generate data from GTF/FASTA sources."
                ) from e
            raise

    files["genes_db"] = _download_or_raise(genes_db_url, base_path, is_gz=False)
    files["fasta_file"] = _download_or_raise(fasta_gz_url, base_path, is_gz=True)
    files["cons_file"] = None
    for ext in (".db", ".pkl"):
        try:
            files["cons_file"] = download(f"{base_url}/{organism}/conservation{ext}", base_path, skip_existing=skip_existing)
            break
        except Exception:
            continue  # Conservation optional for prebuilt path
    return files


def _shell_config_path() -> Optional[Path]:
    """Path to the user's shell config file (~/.zshrc or ~/.bashrc). Prefer zsh if SHELL contains zsh."""
    home = Path.home()
    shell = os.environ.get("SHELL", "")
    if "zsh" in shell:
        rc = home / ".zshrc"
    else:
        rc = home / ".bashrc"
    if rc.exists():
        return rc
    # If preferred rc doesn't exist, try the other
    other = home / ".zshrc" if rc.name == ".bashrc" else home / ".bashrc"
    return other if other.exists() else rc


def _append_seqmat_env_to_shell(data_root: Path, organism: str) -> None:
    """Append SEQMAT_DATA_DIR and SEQMAT_DEFAULT_ORGANISM to the user's shell config if not already set."""
    rc = _shell_config_path()
    if not rc or not rc.is_file():
        return
    try:
        text = rc.read_text()
    except OSError:
        return
    if "SEQMAT_DATA_DIR" in text:
        return
    block = (
        "\n# SeqMat config-less data root (added by seqmat setup)\n"
        f"export SEQMAT_DATA_DIR=\"{data_root}\"\n"
        f"export SEQMAT_DEFAULT_ORGANISM=\"{organism}\"\n"
    )
    try:
        with open(rc, "a") as f:
            f.write(block)
        print(f"Added SEQMAT_DATA_DIR and SEQMAT_DEFAULT_ORGANISM to {rc}. Restart your shell or run: source {rc}")
    except OSError:
        pass


def setup_genomics_data(
    basepath: str,
    organism: Optional[str] = None,
    force: bool = False,
    pickup: bool = False,
    n_jobs: Optional[int] = None,
    from_prebuilt: bool = True,
) -> None:
    """
    Set up genomics data for a specific organism.

    By default downloads prebuilt genes.db and FASTA from the SeqMat S3 bucket.
    Use from_prebuilt=False (or CLI --build-from-sources) to generate data from
    GTF/FASTA sources instead.

    Args:
        basepath: Base directory for storing genomic data
        organism: Organism identifier ('hg38' or 'mm39')
        force: Force overwrite existing data
        pickup: Resume interrupted setup, reuse existing downloaded files
        n_jobs: Number of parallel workers for GTF parsing (default: CPU count - 1). Use 1 for serial.
        from_prebuilt: If True (default), download prebuilt genes.db and FASTA from S3.
            If False, download GTF/FASTA from Ensembl/UCSC and build genes.db locally.
    """
    if organism is None:
        organism = get_default_organism()
    base_path = Path(basepath) / organism

    # Define paths: single base path (no chromosomes/ or annotations/ subdirs unless needed)
    config_paths = {
        'CHROM_SOURCE': str(base_path),
        'MRNA_PATH': str(base_path),
        'MISSPLICING_PATH': str(base_path / 'missplicing'),
        'ONCOSPLICE_PATH': str(base_path / 'oncosplice'),
        'BASE': str(base_path),
        'TEMP': str(base_path / 'temp')
    }

    # Load existing config
    config = load_config()

    # Check if organism already configured
    if organism in config and not force and not pickup:
        print(f"Organism {organism} already configured. Use force=True to overwrite or pickup=True to resume.")
        return

    # Check if directory exists
    if base_path.exists() and any(base_path.iterdir()) and not force and not pickup:
        print(f"Directory {base_path} not empty. Use force=True to overwrite or pickup=True to resume.")
        return

    # Create directory structure
    print(f"Setting up genomics data in {base_path}")
    base_path.mkdir(parents=True, exist_ok=True)

    if from_prebuilt:
        # Default: download prebuilt genes.db and FASTA from S3
        print(f"Downloading prebuilt data for {organism} from S3...")
        try:
            files = download_prebuilt_data(organism, base_path, skip_existing=pickup)
        except PrebuiltDataUnavailableError as e:
            print(str(e))
            raise
        config_paths['fasta_full_genome'] = str(files['fasta_file'])
        config_paths['genes_db'] = str(files['genes_db'])
    else:
        # Build from sources: download GTF/FASTA, parse GTF → genes.db
        print(f"Downloading source data and building genes.db for {organism}...")
        files = download_genome_data(organism, base_path, skip_existing=pickup)
        config_paths['fasta_full_genome'] = str(files['fasta_file'])
        cons_data = None
        if files['cons_file'] is not None:
            cons_data = load_conservation(files['cons_file'])
        annotation_jobs = n_jobs if n_jobs is not None else max(1, cpu_count() - 1)
        retrieve_and_parse_ensembl_annotations(
            base_path,
            files['ensembl_file'],
            cons_data,
            gtex_file=files['gtex_file'],
            n_jobs=annotation_jobs,
        )
        config_paths['genes_db'] = str(base_path / 'genes.db')
        if not pickup:
            for key, file_path in files.items():
                if key == 'fasta_file':
                    continue
                if file_path and file_path.exists():
                    file_path.unlink()

    # Create additional directories
    for path_key in ['MISSPLICING_PATH', 'ONCOSPLICE_PATH', 'TEMP']:
        Path(config_paths[path_key]).mkdir(parents=True, exist_ok=True)

    # Save config file only when not using config-less mode (SEQMAT_DATA_DIR)
    from .config import get_data_base, CONFIG_FILE
    if get_data_base() is None:
        config[organism] = config_paths
        save_config(config)
        print(f"Configuration saved to: {CONFIG_FILE}")
    else:
        print("SEQMAT_DATA_DIR is set; no config file written. Point it at this path to use the data.")
    _append_seqmat_env_to_shell(Path(basepath).resolve(), organism)
    print(f"Successfully set up genomics data for {organism} in {basepath}")
    print("You can now use Gene.from_file() to load gene data.")


def set_fasta_path(fasta_path: str, organism: Optional[str] = None) -> None:
    """
    Set the full genome FASTA path for an organism.

    This is useful when you have a FASTA file but haven't run the full setup,
    or need to update the FASTA path.

    Args:
        fasta_path: Path to the full genome FASTA file
        organism: Organism identifier (uses default if None)
    """
    if organism is None:
        organism = get_default_organism()

    fasta_path = Path(fasta_path)
    if not fasta_path.exists():
        raise ValueError(f"FASTA file not found: {fasta_path}")

    config = load_config()

    if organism not in config:
        config[organism] = {}

    config[organism]['fasta_full_genome'] = str(fasta_path)
    save_config(config)

    from .config import CONFIG_FILE
    print(f"Set fasta_full_genome for {organism} to: {fasta_path}")
    print(f"Configuration saved to: {CONFIG_FILE}")


def list_available_organisms() -> List[str]:
    """
    List all organisms that have been configured.
    
    Returns:
        List of organism names (e.g., ['hg38', 'mm39'])
    """
    config = load_config()
    # Filter out non-organism keys
    non_organism_keys = {'default_organism', 'directory_structure', 'organism_data'}
    organisms = [k for k in config.keys() if k not in non_organism_keys and isinstance(config[k], dict)]
    return organisms


def list_supported_organisms() -> List[str]:
    """
    List all organisms that can be downloaded and configured.
    
    Returns:
        List of supported organism names
    """
    return get_available_organisms()


def get_organism_info(organism: str) -> Dict[str, Any]:
    """
    Get detailed information about a configured organism.
    
    Args:
        organism: Organism identifier (e.g., 'hg38')
        
    Returns:
        Dictionary with organism information
    """
    try:
        from .config import get_organism_config
        config = get_organism_config(organism)
    except ValueError:
        return {"error": f"Organism '{organism}' not configured"}
    
    info = {
        "organism": organism,
        "configured": True,
        "paths": {k: str(v) for k, v in config.items()},
        "data_available": {}
    }
    
    # Check what data is actually available
    mrna_path = config.get('MRNA_PATH')
    if mrna_path and Path(mrna_path).exists():
        biotype_dirs = [d for d in Path(mrna_path).iterdir() if d.is_dir()]
        info["data_available"]["biotypes"] = [d.name for d in biotype_dirs]
        
        # Count genes per biotype
        gene_counts = {}
        for biotype_dir in biotype_dirs:
            gene_files = list(biotype_dir.glob("*.pkl"))
            gene_counts[biotype_dir.name] = len(gene_files)
        info["data_available"]["gene_counts"] = gene_counts
    
    # Check chromosome data
    chrom_path = config.get('CHROM_SOURCE')
    if chrom_path and Path(chrom_path).exists():
        chrom_files = list(Path(chrom_path).glob("*.fasta"))
        info["data_available"]["chromosomes"] = [f.stem for f in chrom_files]
    
    return info


def list_gene_biotypes(organism: str) -> List[str]:
    """
    List all gene biotypes available for an organism.
    
    Args:
        organism: Organism identifier (e.g., 'hg38')
        
    Returns:
        List of biotype names (e.g., ['protein_coding', 'lncRNA'])
    """
    try:
        from .config import get_organism_config
        config = get_organism_config(organism)
    except ValueError:
        return []
    
    mrna_path = config.get('MRNA_PATH')
    if not mrna_path or not Path(mrna_path).exists():
        return []
    
    biotype_dirs = [d.name for d in Path(mrna_path).iterdir() if d.is_dir()]
    return sorted(biotype_dirs)


def count_genes(organism: str, biotype: Optional[str] = None) -> Dict[str, int]:
    """
    Count genes for an organism, optionally filtered by biotype.
    
    Args:
        organism: Organism identifier (e.g., 'hg38')
        biotype: Optional biotype filter (e.g., 'protein_coding')
        
    Returns:
        Dictionary with gene counts per biotype
    """
    try:
        from .config import get_organism_config
        config = get_organism_config(organism)
    except ValueError:
        return {}
    
    mrna_path = config.get('MRNA_PATH')
    if not mrna_path or not Path(mrna_path).exists():
        return {}
    
    gene_counts = {}
    
    if biotype:
        # Count genes for specific biotype
        biotype_path = Path(mrna_path) / biotype
        if biotype_path.exists():
            gene_files = list(biotype_path.glob("*.pkl"))
            gene_counts[biotype] = len(gene_files)
    else:
        # Count genes for all biotypes
        for biotype_dir in Path(mrna_path).iterdir():
            if biotype_dir.is_dir():
                gene_files = list(biotype_dir.glob("*.pkl"))
                gene_counts[biotype_dir.name] = len(gene_files)
    
    return gene_counts


def get_gene_list(organism: str, biotype: str, limit: Optional[int] = None) -> List[str]:
    """
    Get list of gene names for a specific organism and biotype.
    
    Args:
        organism: Organism identifier (e.g., 'hg38')
        biotype: Gene biotype (e.g., 'protein_coding')
        limit: Optional limit on number of genes to return
        
    Returns:
        List of gene names
    """
    try:
        from .config import get_organism_config
        config = get_organism_config(organism)
    except ValueError:
        return []
    
    mrna_path = config.get('MRNA_PATH')
    if not mrna_path:
        return []
    
    biotype_path = Path(mrna_path) / biotype
    if not biotype_path.exists():
        return []
    
    gene_files = list(biotype_path.glob("*.pkl"))
    gene_names = []
    
    for gene_file in gene_files:
        # Extract gene name from filename: mrnas_ENSG123_GENENAME.pkl
        filename = gene_file.stem
        parts = filename.split('_')
        if len(parts) >= 3:
            gene_name = '_'.join(parts[2:])  # Handle gene names with underscores
            gene_names.append(gene_name)
    
    gene_names.sort()
    
    if limit:
        gene_names = gene_names[:limit]
    
    return gene_names


def data_summary() -> Dict[str, Any]:
    """
    Get a comprehensive summary of all available genomics data.
    
    Returns:
        Dictionary with complete data summary
    """
    summary = {
        "supported_organisms": list_supported_organisms(),
        "configured_organisms": list_available_organisms(),
        "organisms": {}
    }
    
    # Get detailed info for each configured organism
    for organism in summary["configured_organisms"]:
        try:
            summary["organisms"][organism] = get_organism_info(organism)
        except Exception as e:
            summary["organisms"][organism] = {"error": f"Configuration error: {str(e)}"}
    
    # Add total counts
    total_genes = 0
    total_biotypes = set()
    
    for org_info in summary["organisms"].values():
        if "gene_counts" in org_info.get("data_available", {}):
            for biotype, count in org_info["data_available"]["gene_counts"].items():
                total_genes += count
                total_biotypes.add(biotype)
    
    summary["totals"] = {
        "organisms": len(summary["configured_organisms"]),
        "biotypes": len(total_biotypes),
        "genes": total_genes
    }
    
    return summary


def print_data_summary():
    """Print a formatted summary of all available genomics data."""
    summary = data_summary()
    
    print("🧬 SeqMat Genomics Data Summary")
    print("=" * 40)
    
    # Totals
    totals = summary["totals"]
    print(f"📊 Total: {totals['organisms']} organisms, {totals['biotypes']} biotypes, {totals['genes']} genes")
    print()
    
    # Supported organisms
    print("🌍 Supported Organisms:")
    supported = list_supported_organisms()
    for org in supported:
        status = "✅ Configured" if org in summary["configured_organisms"] else "❌ Not configured"
        name = DEFAULT_ORGANISM_DATA.get(org, {}).get('name', org)
        print(f"   {org}: {name} - {status}")
    print()
    
    # Detailed organism info
    for organism, info in summary["organisms"].items():
        print(f"📁 {organism.upper()} Data:")
        
        if "error" in info:
            print(f"   ❌ {info['error']}")
            continue
        
        data_avail = info.get("data_available", {})
        
        # Gene counts by biotype
        if "gene_counts" in data_avail:
            print("   Gene Types:")
            for biotype, count in sorted(data_avail["gene_counts"].items()):
                print(f"     {biotype}: {count:,} genes")
        
        # Chromosomes
        if "chromosomes" in data_avail:
            chroms = data_avail["chromosomes"]
            print(f"   Chromosomes: {len(chroms)} available ({', '.join(sorted(chroms)[:5])}{'...' if len(chroms) > 5 else ''})")
        
        # Paths
        print("   Data Paths:")
        for path_name, path_value in info["paths"].items():
            exists = "✅" if Path(path_value).exists() else "❌"
            print(f"     {path_name}: {exists} {path_value}")
        
        print()


def get_all_genes(organism: str, biotype: Optional[str] = None) -> List[Dict[str, str]]:
    """
    Get all available genes for an organism, optionally filtered by biotype.
    
    Args:
        organism: Organism identifier (e.g., 'hg38')
        biotype: Optional biotype filter (e.g., 'protein_coding', 'lncRNA')
        
    Returns:
        List of dictionaries with gene information including:
        - organism: The organism identifier
        - biotype: The gene biotype
        - gene_name: The gene name
        - gene_id: The gene ID (e.g., ENSG00000133703)
    """
    results = []
    
    try:
        from .config import get_organism_config
        config = get_organism_config(organism)
    except ValueError:
        return []
    
    mrna_path = config.get('MRNA_PATH')
    if not mrna_path or not Path(mrna_path).exists():
        return []
    
    # Get biotypes to search
    if biotype:
        biotypes = [biotype]
    else:
        biotypes = list_gene_biotypes(organism)
    
    # Collect all genes from each biotype
    for bt in biotypes:
        biotype_path = Path(mrna_path) / bt
        if not biotype_path.exists():
            continue
        
        gene_files = list(biotype_path.glob("*.pkl"))
        
        for gene_file in gene_files:
            # Extract gene info from filename: mrnas_ENSG123_GENENAME.pkl
            filename = gene_file.stem
            parts = filename.split('_', 2)  # Split into at most 3 parts
            if len(parts) >= 3:
                gene_id = parts[1]
                gene_name = parts[2]
                
                results.append({
                    "organism": organism,
                    "biotype": bt,
                    "gene_name": gene_name,
                    "gene_id": gene_id
                })
    
    # Sort results by gene name for consistency
    results.sort(key=lambda x: x['gene_name'])
    
    return results


def search_genes(organism: str, query: str, biotype: Optional[str] = None, limit: int = 10) -> List[Dict[str, str]]:
    """
    Search for genes by name or ID pattern.
    
    Args:
        organism: Organism identifier (e.g., 'hg38')
        query: Search query (gene name or ID pattern)
        biotype: Optional biotype filter
        limit: Maximum number of results
        
    Returns:
        List of dictionaries with gene information
    """
    results = []
    query_upper = query.upper()  # Gene names are typically uppercase
    query_lower = query.lower()
    
    try:
        from .config import get_organism_config
        config = get_organism_config(organism)
    except ValueError:
        return []
    
    mrna_path = config.get('MRNA_PATH')
    if not mrna_path:
        return []
    
    if biotype:
        biotypes = [biotype]
    else:
        biotypes = list_gene_biotypes(organism)
    
    for bt in biotypes:
        biotype_path = Path(mrna_path) / bt
        if not biotype_path.exists():
            continue
        
        gene_files = list(biotype_path.glob("*.pkl"))
        
        for gene_file in gene_files:
            # Extract gene info from filename: mrnas_ENSG123_GENENAME.pkl
            filename = gene_file.stem
            parts = filename.split('_', 2)  # Split into at most 3 parts
            if len(parts) >= 3:
                gene_id = parts[1]
                gene_name = parts[2]
                
                # Search in both gene name and gene ID
                if (query_upper in gene_name or 
                    query_lower in gene_name.lower() or
                    query in gene_id or
                    query_lower in gene_id.lower()):
                    results.append({
                        "organism": organism,
                        "biotype": bt,
                        "gene_name": gene_name,
                        "gene_id": gene_id
                    })
                    
                    if len(results) >= limit:
                        return results
    
    return results


def test_installation(organism: str = None, verbose: bool = True) -> Dict[str, Any]:
    """
    Comprehensive test of SeqMat installation and data setup.
    
    Args:
        organism: Organism to test (uses default if None)
        verbose: Print detailed output
        
    Returns:
        Dictionary with test results
    """
    from .seqmat import SeqMat
    from .gene import Gene
    from .config import get_organism_config
    
    if organism is None:
        organism = get_default_organism()
    
    results = {
        "organism": organism,
        "tests_passed": 0,
        "tests_failed": 0,
        "errors": [],
        "warnings": []
    }
    
    def test_step(name: str, func):
        """Execute a test step and record results"""
        if verbose:
            print(f"\n🧪 Testing {name}...", end="", flush=True)
        try:
            result = func()
            results["tests_passed"] += 1
            if verbose:
                print(" ✅")
            return True, result
        except Exception as e:
            results["tests_failed"] += 1
            results["errors"].append(f"{name}: {str(e)}")
            if verbose:
                print(f" ❌ {str(e)}")
            return False, None
    
    if verbose:
        print(f"🔬 Running comprehensive tests for {organism}")
        print("=" * 50)
    
    # Test 1: Configuration
    def test_config():
        config = get_organism_config(organism)
        assert config, "No configuration found"
        assert 'MRNA_PATH' in config, "Missing MRNA_PATH"
        assert 'CHROM_SOURCE' in config, "Missing CHROM_SOURCE"
        return config
    
    success, config = test_step("Configuration", test_config)
    if not success:
        return results
    
    # Test 2: Directory structure
    def test_directories():
        paths_exist = {}
        for key, path in config.items():
            exists = Path(path).exists()
            paths_exist[key] = exists
            if not exists and key in ['MRNA_PATH', 'CHROM_SOURCE']:
                raise ValueError(f"Critical path missing: {key} -> {path}")
        return paths_exist
    
    test_step("Directory structure", test_directories)
    
    # Test 3: Chromosome files
    def test_chromosomes():
        chrom_path = Path(config['CHROM_SOURCE'])
        chrom_files = list(chrom_path.glob("*.fasta"))
        if not chrom_files:
            raise ValueError("No chromosome files found")
        # Test we can read from first chromosome
        first_chrom = chrom_files[0].stem
        return {"count": len(chrom_files), "first": first_chrom}
    
    success, chrom_info = test_step("Chromosome files", test_chromosomes)
    
    # Test 4: Gene annotations
    def test_annotations():
        biotypes = list_gene_biotypes(organism)
        if not biotypes:
            raise ValueError("No gene biotypes found")
        gene_counts = count_genes(organism)
        total_genes = sum(gene_counts.values())
        if total_genes == 0:
            raise ValueError("No genes found")
        return {"biotypes": len(biotypes), "total_genes": total_genes}
    
    success, gene_info = test_step("Gene annotations", test_annotations)
    
    # Test 5: Gene search
    def test_gene_search():
        # Search for common genes
        test_genes = ['KRAS', 'TP53', 'EGFR', 'BRCA1', 'MYC']
        found_genes = []
        for gene in test_genes:
            results = search_genes(organism, gene, limit=1)
            if results:
                found_genes.append(gene)
        if not found_genes:
            raise ValueError("Could not find any common genes")
        return {"searched": len(test_genes), "found": len(found_genes)}
    
    test_step("Gene search", test_gene_search)
    
    # Test 6: Load a specific gene
    def test_gene_loading():
        # Try to find and load a protein-coding gene
        protein_genes = get_gene_list(organism, 'protein_coding', limit=10)
        if not protein_genes:
            raise ValueError("No protein coding genes found")
        
        # Try to load the first gene
        gene_name = protein_genes[0]
        gene = Gene.from_file(gene_name, organism=organism)
        
        # Verify gene properties
        assert gene.gene_name == gene_name, f"Gene name mismatch: {gene.gene_name} != {gene_name}"
        assert len(gene.transcripts) > 0, "Gene has no transcripts"
        
        return {
            "gene_name": gene_name,
            "gene_id": gene.gene_id,
            "transcripts": len(gene.transcripts),
            "chromosome": gene.chrm
        }
    
    success, gene_data = test_step("Gene loading", test_gene_loading)
    
    # Test 7: SeqMat FASTA loading
    def test_fasta_loading():
        if not chrom_info:
            raise ValueError("No chromosome info available")
        
        # Use first available chromosome
        chrom = chrom_info['first']
        # Try to load a small region
        seq = SeqMat.from_fasta(
            genome=organism,
            chrom=chrom,
            start=1000000,
            end=1000100
        )
        
        assert len(seq) == 101, f"Expected 101bp, got {len(seq)}"
        assert seq.seq, "No sequence loaded"
        assert all(c in 'ACGTN' for c in seq.seq), "Invalid nucleotides in sequence"
        
        return {
            "chromosome": chrom,
            "length": len(seq),
            "sequence_preview": seq.seq[:20] + "..."
        }
    
    if chrom_info:
        test_step("FASTA loading", test_fasta_loading)
    
    # Test 8: Transcript operations
    def test_transcript_ops():
        if not gene_data:
            raise ValueError("No gene data available")
        
        # Reload the gene
        gene = Gene.from_file(gene_data['gene_name'], organism=organism)
        
        # Get primary transcript
        transcript = gene.transcript()
        
        # Test basic properties
        assert transcript.transcript_id, "No transcript ID"
        assert len(transcript.exons) > 0, "No exons found"
        
        # Try to generate mature mRNA
        transcript.generate_mature_mrna()
        assert transcript.mature_mrna is not None, "Failed to generate mature mRNA"
        assert len(transcript.mature_mrna) > 0, "Empty mature mRNA"
        
        return {
            "transcript_id": transcript.transcript_id,
            "exons": len(transcript.exons),
            "mature_mrna_length": len(transcript.mature_mrna) if transcript.mature_mrna else 0,
            "protein_coding": transcript.protein_coding
        }
    
    if gene_data:
        test_step("Transcript operations", test_transcript_ops)
    
    # Test 9: SeqMat mutations and complement
    def test_mutations():
        # Create a simple sequence
        seq = SeqMat("ATCGATCGATCG")
        original_seq = seq.seq
        
        # Test complement
        comp = seq.complement()
        expected_comp = "TAGCTAGCTAGC"
        assert comp.seq == expected_comp, f"Complement failed: got {comp.seq}, expected {expected_comp}"
        
        # Test reverse complement
        rev_comp = seq.clone()
        rev_comp.reverse_complement()
        expected_rev_comp = "CGATCGATCGAT"
        assert rev_comp.seq == expected_rev_comp, f"Reverse complement failed: got {rev_comp.seq}, expected {expected_rev_comp}"
        
        # Apply mutations
        mutations = [
            (3, "C", "G"),      # SNP
            (6, "-", "AAA"),    # Insertion
            (10, "TC", "-")     # Deletion
        ]
        seq.apply_mutations(mutations)
        
        assert len(seq.mutations) == 3, f"Expected 3 mutations, got {len(seq.mutations)}"
        assert seq.seq != original_seq, "Sequence unchanged after mutations"
        
        return {
            "original_length": 12,
            "mutated_length": len(seq),
            "mutations_applied": len(seq.mutations),
            "complement_test": "passed",
            "reverse_complement_test": "passed"
        }
    
    test_step("Sequence mutations and complement", test_mutations)
    
    # Test 10: Data summary
    def test_data_summary():
        summary = data_summary()
        assert 'configured_organisms' in summary, "Missing configured organisms"
        assert organism in summary['configured_organisms'], f"{organism} not in configured organisms"
        return {
            "configured_organisms": len(summary['configured_organisms']),
            "total_genes": summary['totals']['genes']
        }
    
    test_step("Data summary", test_data_summary)
    
    # Summary
    if verbose:
        print("\n" + "=" * 50)
        print(f"📊 Test Summary for {organism}:")
        print(f"   ✅ Passed: {results['tests_passed']}")
        print(f"   ❌ Failed: {results['tests_failed']}")
        
        if results['errors']:
            print("\n❌ Errors:")
            for error in results['errors']:
                print(f"   - {error}")
        
        if results['warnings']:
            print("\n⚠️  Warnings:")
            for warning in results['warnings']:
                print(f"   - {warning}")
        
        if results['tests_failed'] == 0:
            print("\n✨ All tests passed! SeqMat is properly installed and configured.")
        else:
            print("\n⚠️  Some tests failed. Please check the errors above.")
    
    return results


def generate_random_sequence(k=15_000):
    import random
    bases = ['A', 'C', 'G', 'T']
    seq = ''.join(random.choices(bases, k=k))
    return seq


def available_genes(organism: str = 'hg38'):
    """
    Yield gene names (symbols) found in the protein_coding directory for an organism.

    This is useful for quickly checking which genes have data files available
    that can be loaded with Gene.from_file().

    Args:
        organism: Organism identifier (e.g., 'hg38', 'mm39')

    Yields:
        Gene symbols (e.g., 'KRAS', 'TP53', 'BRCA1')

    Example:
        >>> genes = list(available_genes('hg38'))
        >>> 'KRAS' in genes
        True
        >>> for gene in available_genes():
        ...     print(gene)
    """
    import os
    from .config import get_organism_config

    config = get_organism_config(organism)
    mrna_path = config.get('MRNA_PATH')

    if not mrna_path:
        raise ValueError(f"MRNA_PATH not found in config for organism '{organism}'")

    protein_coding_path = Path(mrna_path) / 'protein_coding'

    if not protein_coding_path.exists():
        raise ValueError(f"protein_coding directory not found at {protein_coding_path}")

    for file in os.listdir(protein_coding_path):
        if file.endswith('.pkl'):
            gene = file.split('_')[-1].removesuffix('.pkl')
            yield gene