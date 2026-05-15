"""HTTP downloads for genome data: GTF, FASTA, and prebuilt assets.

Public:
- :func:`download`, :func:`download_and_ungzip`
- :func:`download_genome_data`, :func:`download_prebuilt_data`
- :class:`PrebuiltDataUnavailableError` (alias of :class:`DataUnavailableError`)
"""
from __future__ import annotations

import gzip
import logging
import shutil
from pathlib import Path
from typing import Dict, Optional

import requests
from tqdm import tqdm

from .config import DEFAULT_ORGANISM_DATA, get_prebuilt_data_base_url
from .errors import DataUnavailableError, PrebuiltDataUnavailableError

_log = logging.getLogger(__name__)


def download(external_url: str, local_path: Path, skip_existing: bool = False) -> Path:
    """Download a file from URL to local path. Streams in 8 KiB chunks with a progress bar."""
    local_file = Path(external_url).name
    local_file_path = Path(local_path) / local_file

    if skip_existing and local_file_path.exists():
        _log.info("File already exists: %s", local_file_path.name)
        return local_file_path

    _log.info("Downloading %s", external_url)
    response = requests.get(external_url, stream=True)
    response.raise_for_status()
    total_size = int(response.headers.get("content-length", 0))
    block_size = 8192
    with open(local_file_path, "wb") as f, tqdm(total=total_size, unit="iB", unit_scale=True) as pbar:
        for chunk in response.iter_content(block_size):
            if chunk:
                f.write(chunk)
                pbar.update(len(chunk))
    return local_file_path


def download_and_ungzip(external_url: str, local_path: Path, skip_existing: bool = False) -> Path:
    """Download a gzipped file and decompress it, leaving only the plain file behind."""
    local_file = Path(external_url).name
    output_file = Path(local_path) / local_file.rstrip(".gz")

    if skip_existing and output_file.exists():
        _log.info("Uncompressed file already exists: %s", output_file.name)
        return output_file

    gzipped_file = download(external_url, local_path, skip_existing=skip_existing)
    if gzipped_file.suffix == ".gz" and not output_file.exists():
        _log.info("Decompressing %s", gzipped_file.name)
        with gzip.open(gzipped_file, "rb") as f_in, open(output_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        gzipped_file.unlink()
    return output_file


def _validate_prebuilt_file(
    path: Path,
    min_size: int,
    kind: str,
    url: str,
    allow_fasta_header: bool = False,
) -> None:
    """Raise :class:`DataUnavailableError` and remove ``path`` if it looks like an error page."""
    size = path.stat().st_size
    if size < min_size:
        path.unlink()
        raise DataUnavailableError(
            f"Prebuilt {kind} from {url} is too small ({size} bytes). "
            "The server may have returned an error page. Check the bucket URL, network, and permissions. "
            "Alternatively run setup with --build-from-sources."
        )
    with open(path, "rb") as f:
        head = f.read(64)
    if head.startswith(b"<?") or head.startswith(b"<!") or b"<html" in head[:32].lower():
        path.unlink()
        raise DataUnavailableError(
            f"Prebuilt {kind} from {url} looks like an HTML/XML error page, not data. "
            "Check the bucket URL and permissions. Alternatively run setup with --build-from-sources."
        )
    if allow_fasta_header and not head.lstrip().startswith(b">"):
        path.unlink()
        raise DataUnavailableError(
            f"Prebuilt {kind} from {url} does not look like a FASTA file. "
            "Check the bucket. Alternatively run setup with --build-from-sources."
        )


def download_genome_data(organism: str, base_path: Path, skip_existing: bool = False) -> Dict[str, Optional[Path]]:
    """Download GTF (Ensembl) + FASTA (UCSC) + optional GTEx/conservation for an organism."""
    if organism not in DEFAULT_ORGANISM_DATA:
        raise ValueError(
            f"Organism {organism} not supported for download. Available: {list(DEFAULT_ORGANISM_DATA.keys())}"
        )

    urls = DEFAULT_ORGANISM_DATA[organism].get("urls", {})
    if not urls:
        raise ValueError(f"No URLs configured for organism {organism}")

    url_mapping = {"fasta": "fasta_url", "gtf": "ensembl_url", "gtex": "expression_url"}
    mapped = {url_mapping[k]: v for k, v in urls.items() if k in url_mapping}
    if not mapped:
        raise ValueError(f"No valid URLs found for organism {organism}. Available keys: {list(urls.keys())}")

    files: Dict[str, Optional[Path]] = {}
    if not mapped.get("ensembl_url"):
        raise ValueError(f"No GTF URL configured for organism {organism}")
    files["ensembl_file"] = download_and_ungzip(mapped["ensembl_url"], base_path, skip_existing=skip_existing)

    files["cons_file"] = None
    base_url = get_prebuilt_data_base_url()
    for ext in (".db", ".pkl"):
        try:
            files["cons_file"] = download(
                f"{base_url}/{organism}/conservation{ext}", base_path, skip_existing=skip_existing
            )
            break
        except Exception:
            continue

    files["gtex_file"] = (
        download_and_ungzip(mapped["expression_url"], base_path, skip_existing=skip_existing)
        if mapped.get("expression_url")
        else None
    )

    if not mapped.get("fasta_url"):
        raise ValueError(f"No FASTA URL configured for organism {organism}")
    files["fasta_file"] = download_and_ungzip(mapped["fasta_url"], base_path, skip_existing=skip_existing)
    return files


def download_prebuilt_data(
    organism: str,
    base_path: Path,
    skip_existing: bool = False,
) -> Dict[str, Optional[Path]]:
    """Download prebuilt ``genes.db`` and FASTA from the SeqMat S3 bucket.

    Layout in the bucket: ``{base}/{organism}/genes.db`` and ``{base}/{organism}/{organism}.fa.gz``.
    Conservation is not downloaded here; it is only fetched during ``--build-from-sources``.
    """
    if organism not in DEFAULT_ORGANISM_DATA:
        raise ValueError(
            f"Organism {organism} not supported. Available: {list(DEFAULT_ORGANISM_DATA.keys())}"
        )
    base_path = Path(base_path)
    base_path.mkdir(parents=True, exist_ok=True)
    base_url = get_prebuilt_data_base_url()
    genes_db_url = f"{base_url}/{organism}/genes.db"
    fasta_gz_url = f"{base_url}/{organism}/{organism}.fa.gz"
    fasta_plain_url = f"{base_url}/{organism}/{organism}.fa"
    files: Dict[str, Optional[Path]] = {}

    def _download_or_raise(url: str, dest_dir: Path, is_gz: bool = False) -> Path:
        try:
            if is_gz:
                return download_and_ungzip(url, dest_dir, skip_existing=skip_existing)
            return download(url, dest_dir, skip_existing=skip_existing)
        except requests.HTTPError as e:
            if e.response is not None and e.response.status_code == 404:
                raise DataUnavailableError(
                    f"Prebuilt data not found at {url}. "
                    "Upload genes.db and FASTA to the S3 bucket, or run setup with "
                    "--build-from-sources to generate data from GTF/FASTA sources."
                ) from e
            raise
        except Exception as e:
            if "404" in str(e).lower() or "not found" in str(e).lower():
                raise DataUnavailableError(
                    f"Prebuilt data not found for {organism}. "
                    "Run setup with --build-from-sources to generate data from GTF/FASTA sources."
                ) from e
            raise

    files["genes_db"] = _download_or_raise(genes_db_url, base_path, is_gz=False)
    _validate_prebuilt_file(files["genes_db"], min_size=100_000, kind="genes.db", url=genes_db_url)

    fasta_plain_path = base_path / f"{organism}.fa"
    if skip_existing and fasta_plain_path.exists():
        _log.info("Uncompressed FASTA already exists: %s", fasta_plain_path.name)
        files["fasta_file"] = fasta_plain_path
    else:
        try:
            files["fasta_file"] = _download_or_raise(fasta_gz_url, base_path, is_gz=True)
        except (OSError, requests.HTTPError) as e:
            is_404_403 = (
                isinstance(e, requests.HTTPError)
                and e.response is not None
                and e.response.status_code in (403, 404)
            )
            is_not_gzip = isinstance(e, OSError) and (
                "not a gzipped" in str(e).lower() or "bad gzip" in str(e).lower()
            )
            if not (is_404_403 or is_not_gzip):
                raise
            gz_path = base_path / f"{organism}.fa.gz"
            if gz_path.exists():
                gz_path.unlink()
            files["fasta_file"] = _download_or_raise(fasta_plain_url, base_path, is_gz=False)

    _validate_prebuilt_file(
        files["fasta_file"],
        min_size=1_000_000,
        kind="FASTA",
        url=f"{base_url}/{organism}/",
        allow_fasta_header=True,
    )
    files["cons_file"] = None
    return files


__all__ = [
    "download",
    "download_and_ungzip",
    "download_genome_data",
    "download_prebuilt_data",
    "PrebuiltDataUnavailableError",
    "DataUnavailableError",
]
