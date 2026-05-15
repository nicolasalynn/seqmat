# SeqMat — Data Setup and Configuration

Most users just run `seqmat setup` once and never look at this file. The notes below cover the rest: custom organisms, mirrors, ephemeral environments, multi-user installs, and the full configuration system.

## TL;DR

```bash
pip install seqmat
seqmat setup                              # downloads prebuilt hg38 data (~4 GB)
```

`seqmat setup` appends two lines to your shell rc (`SEQMAT_DATA_DIR`, `SEQMAT_DEFAULT_ORGANISM`). After that, future sessions are config-less.

## Configuration resolution

SeqMat resolves data paths in this order:

1. **Config-less mode (recommended).** If `SEQMAT_DATA_DIR` is set, all paths are derived: `{SEQMAT_DATA_DIR}/{organism}/genes.db`, `{organism}.fa`, etc. `SEQMAT_DEFAULT_ORGANISM` picks the default. No config file is read.
2. **Config file mode.** When `SEQMAT_DATA_DIR` is not set, paths come from `config.json`. Its location is resolved as:
   - `SEQMAT_CONFIG_FILE` → exact path
   - `SEQMAT_CONFIG_DIR` → directory; config is `{dir}/config.json`
   - default: platform user config dir (`~/.config/seqmat/config.json` on Linux)

Inspect what's active:

```python
from seqmat import get_data_base, get_config_file
print(get_data_base())     # Path if config-less mode is active, else None
print(get_config_file())   # Active config.json path
```

## Programmatic setup

```python
from seqmat import setup_genomics_data

# Default: download prebuilt genes.db + FASTA from S3 (no rebuild)
setup_genomics_data("/data/seqmat", organism="hg38")
setup_genomics_data("/data/seqmat", organism="mm39")

# Build genes.db from a fresh GTF (Ensembl 111 by default for hg38)
setup_genomics_data("/data/seqmat", organism="hg38", from_prebuilt=False, n_jobs=8)
```

## CLI reference

```bash
seqmat setup [--path PATH] [--organism hg38|mm39] [--force] [--build-from-sources]
seqmat organisms                                       # status of supported organisms
seqmat summary                                         # data overview
seqmat info     --organism hg38
seqmat biotypes --organism hg38
seqmat count    --organism hg38 [--biotype protein_coding]
seqmat list     --organism hg38 --biotype protein_coding [--limit N]
seqmat search   --organism hg38 --query KRAS [--biotype TYPE] [--limit N]
```

## Directory layout

```
{SEQMAT_DATA_DIR}/
├── hg38/
│   ├── genes.db              # SQLite, gene records keyed by gene_name / gene_id
│   ├── gene_locations.npz    # per-chromosome position index (built lazily on first lookup)
│   ├── hg38.fa               # full genome FASTA
│   ├── hg38.fa.fai           # samtools FASTA index
│   ├── conservation.db       # optional, populated by --build-from-sources
│   ├── missplicing/
│   ├── oncosplice/
│   └── temp/
└── mm39/  (same layout)
```

## Mirroring the prebuilt bucket

By default SeqMat downloads from the public `seqmat-prebuilt-public` S3 bucket. To use a private mirror:

```bash
export SEQMAT_PREBUILT_DATA_BASE_URL="https://your-bucket.s3.eu-north-1.amazonaws.com"
```

Or set `prebuilt_data_base_url` in `config.json`. Expected layout in the bucket:

```
{base}/{organism}/genes.db
{base}/{organism}/{organism}.fa.gz   # or .fa
{base}/{organism}/conservation.db    # optional
```

## Adding a new organism

For organisms beyond hg38/mm39 you have three options.

### 1. Extend the config (no code change)

```python
from seqmat.config import load_config, save_config

config = load_config()
config["dm6"] = {
    "name": "Drosophila melanogaster",
    "urls": {
        "fasta": "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz",
        "gtf":   "https://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf.gz",
    },
}
save_config(config)
setup_genomics_data("/data/seqmat", organism="dm6", from_prebuilt=False)
```

### 2. Add to `DEFAULT_ORGANISM_DATA` in `seqmat/config.py`

Permanent additions (e.g. to upstream a new organism) belong here. Same shape as the dict above.

### 3. Bring-your-own data

Drop a `genes.db` (same schema produced by `--build-from-sources`) and a FASTA into the organism directory, then point config there:

```python
from pathlib import Path
from seqmat.config import load_config, save_config

base = Path("/data/custom_genome")
config = load_config()
config["custom_genome"] = {
    "BASE": str(base),
    "MRNA_PATH": str(base),
    "CHROM_SOURCE": str(base),
    "genes_db": str(base / "genes.db"),
    "fasta_full_genome": str(base / "custom_genome.fa"),
}
save_config(config)
```

## Environment recipes

**Ephemeral runners (Docker, Run.ai, GitHub Actions):**

```bash
export SEQMAT_DATA_DIR=/mnt/persistent/seqmat
export SEQMAT_DEFAULT_ORGANISM=hg38
```

Mount the persistent volume so `genes.db` and the FASTA survive between runs.

**Shared multi-user system:**

```bash
# /etc/profile.d/seqmat.sh
export SEQMAT_DATA_DIR=/shared/genomics/seqmat
export SEQMAT_DEFAULT_ORGANISM=hg38
```

Everyone reads from the same path; nobody re-downloads.

**Custom config location (per-user override):**

```bash
export SEQMAT_CONFIG_FILE=/path/to/my/config.json
# or
export SEQMAT_CONFIG_DIR=/path/to/config/dir   # reads {dir}/config.json
```

## Troubleshooting

| Symptom                                 | Likely cause                                                |
| --------------------------------------- | ----------------------------------------------------------- |
| `Organism 'X' not configured`           | Run `seqmat setup --organism X`, or check `seqmat organisms` |
| `Gene.from_file` returns `None`         | Gene name not in `genes.db` for that organism. Try `search_genes("X", "KRAS")` |
| `FileNotFoundError` on FASTA            | `fasta_full_genome` path wrong; run `set_fasta_path("/path/to/hg38.fa")` |
| Download failures                       | Check network / bucket URL / `SEQMAT_PREBUILT_DATA_BASE_URL` |
| `PrebuiltDataUnavailableError`          | Bucket returned 404; retry with `--build-from-sources`      |
| Permission errors                       | Pick a `SEQMAT_DATA_DIR` you can write to                   |

Disk-space rule of thumb: hg38 ≈ 4 GB, mm39 ≈ 3 GB.

## Configuration API

```python
from seqmat.config import (
    load_config, save_config,
    get_default_organism, get_available_organisms,
    get_organism_config, get_organism_info,
    get_prebuilt_data_base_url,
)
```

These are the building blocks the CLI uses; reach for them when you want to script setup or interrogate state from Python.
