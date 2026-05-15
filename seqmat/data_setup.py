"""High-level ``setup_genomics_data`` entry point and helpers.

The :mod:`download`, :mod:`gtf`, and :mod:`conservation` modules do the work;
this module wires them together for the ``seqmat setup`` CLI command.
"""
from __future__ import annotations

import logging
import os
from multiprocessing import cpu_count
from pathlib import Path
from typing import Optional

from .config import (
    CONFIG_FILE,
    get_data_base,
    get_default_organism,
    load_config,
    save_config,
)
from .conservation import load_conservation
from .download import (
    PrebuiltDataUnavailableError,
    download_genome_data,
    download_prebuilt_data,
)
from .gtf import retrieve_and_parse_ensembl_annotations

_log = logging.getLogger(__name__)


def _shell_config_path() -> Optional[Path]:
    """Path to the user's shell config file. Prefer ``~/.zshrc`` when ``$SHELL`` contains ``zsh``."""
    home = Path.home()
    shell = os.environ.get("SHELL", "")
    rc = home / (".zshrc" if "zsh" in shell else ".bashrc")
    if rc.exists():
        return rc
    other = home / ".zshrc" if rc.name == ".bashrc" else home / ".bashrc"
    return other if other.exists() else rc


def _append_seqmat_env_to_shell(data_root: Path, organism: str) -> None:
    """Append ``SEQMAT_DATA_DIR`` / ``SEQMAT_DEFAULT_ORGANISM`` to the user's shell rc, idempotently."""
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
        f'export SEQMAT_DATA_DIR="{data_root}"\n'
        f'export SEQMAT_DEFAULT_ORGANISM="{organism}"\n'
    )
    try:
        with open(rc, "a") as f:
            f.write(block)
        _log.info(
            "Added SEQMAT_DATA_DIR and SEQMAT_DEFAULT_ORGANISM to %s. "
            "Restart your shell or run: source %s",
            rc, rc,
        )
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
    """Set up genomics data for a specific organism.

    Args:
        basepath: Base directory for storing genomic data.
        organism: Organism identifier (e.g. ``"hg38"`` or ``"mm39"``). Uses the configured
            default if omitted.
        force: Overwrite existing data.
        pickup: Resume an interrupted setup, reusing any already-downloaded files.
        n_jobs: Workers for GTF parsing (default: CPU count - 1). ``1`` is serial.
        from_prebuilt: When ``True`` (default), download prebuilt ``genes.db`` + FASTA from S3.
            When ``False``, download GTF/FASTA from Ensembl/UCSC and build ``genes.db`` locally.
    """
    if organism is None:
        organism = get_default_organism()
    base_path = Path(basepath) / organism

    config_paths = {
        "CHROM_SOURCE": str(base_path),
        "MRNA_PATH": str(base_path),
        "MISSPLICING_PATH": str(base_path / "missplicing"),
        "ONCOSPLICE_PATH": str(base_path / "oncosplice"),
        "BASE": str(base_path),
        "TEMP": str(base_path / "temp"),
    }

    config = load_config()
    if organism in config and not force and not pickup:
        _log.warning(
            "Organism %s already configured. Use force=True to overwrite or pickup=True to resume.",
            organism,
        )
        return

    if base_path.exists() and any(base_path.iterdir()) and not force and not pickup:
        _log.warning(
            "Directory %s not empty. Use force=True to overwrite or pickup=True to resume.",
            base_path,
        )
        return

    _log.info("Setting up genomics data in %s", base_path)
    base_path.mkdir(parents=True, exist_ok=True)

    if from_prebuilt:
        _log.info("Downloading prebuilt data for %s from S3...", organism)
        try:
            files = download_prebuilt_data(organism, base_path, skip_existing=pickup)
        except PrebuiltDataUnavailableError as e:
            _log.error("%s", e)
            raise
        config_paths["fasta_full_genome"] = str(files["fasta_file"])
        config_paths["genes_db"] = str(files["genes_db"])
    else:
        _log.info("Downloading source data and building genes.db for %s...", organism)
        files = download_genome_data(organism, base_path, skip_existing=pickup)
        config_paths["fasta_full_genome"] = str(files["fasta_file"])
        cons_data = None
        if files["cons_file"] is not None:
            try:
                cons_data = load_conservation(files["cons_file"])
            except Exception:
                cons_data = None
        annotation_jobs = n_jobs if n_jobs is not None else max(1, cpu_count() - 1)
        retrieve_and_parse_ensembl_annotations(
            base_path,
            files["ensembl_file"],
            cons_data,
            gtex_file=files["gtex_file"],
            n_jobs=annotation_jobs,
        )
        config_paths["genes_db"] = str(base_path / "genes.db")
        if not pickup:
            for key, file_path in files.items():
                if key == "fasta_file":
                    continue
                if file_path and file_path.exists():
                    file_path.unlink()

    for path_key in ("MISSPLICING_PATH", "ONCOSPLICE_PATH", "TEMP"):
        Path(config_paths[path_key]).mkdir(parents=True, exist_ok=True)

    if get_data_base() is None:
        config[organism] = config_paths
        save_config(config)
        _log.info("Configuration saved to: %s", CONFIG_FILE)
    else:
        _log.info("SEQMAT_DATA_DIR is set; no config file written.")
    _append_seqmat_env_to_shell(Path(basepath).resolve(), organism)
    _log.info("Successfully set up genomics data for %s in %s", organism, basepath)


def set_fasta_path(fasta_path: str, organism: Optional[str] = None) -> None:
    """Set the full genome FASTA path for an organism in the saved config.

    Useful when you have a FASTA but haven't run the full setup, or need to point
    SeqMat at a different reference build.
    """
    if organism is None:
        organism = get_default_organism()
    fa_path = Path(fasta_path)
    if not fa_path.exists():
        raise ValueError(f"FASTA file not found: {fa_path}")

    config = load_config()
    config.setdefault(organism, {})
    config[organism]["fasta_full_genome"] = str(fa_path)
    save_config(config)
    _log.info("Set fasta_full_genome for %s to: %s (config: %s)", organism, fa_path, CONFIG_FILE)


__all__ = ["setup_genomics_data", "set_fasta_path"]
