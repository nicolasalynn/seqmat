# SeqMat

**Lightning-fast genomic sequence matrix library with mutation tracking**

SeqMat is a comprehensive Python library for genomic sequence analysis, providing efficient tools for working with genes, transcripts, and genomic sequences. Features include full mutation tracking, splicing analysis, and high-performance sequence manipulation.

## Key Features

- **Fast Sequence Operations**: Vectorized genomic sequence matrix with efficient slicing and manipulation
- **Comprehensive Mutation Tracking**: SNPs, insertions, deletions with full history and conflict detection
- **Gene & Transcript Analysis**: Load and analyze gene structures, exons, introns, and splice sites
- **Conservation Integration**: Built-in support for conservation scores and sequence analysis
- **Multi-organism Support**: Human (hg38) and mouse (mm39) genome support
- **Data Inspection Tools**: Comprehensive utilities to explore available data, gene counts, and biotypes
- **Command-line Interface**: Full CLI for data management and inspection
- **Expression Data**: Integration with GTEx tissue expression data
- **Flexible Installation**: Core functionality with optional bioinformatics dependencies

## Installation

### Basic Installation
```bash
pip install seqmat
```

After installing, to use genes and transcripts you need to fetch organism data once (see [Data Setup and Configuration](#data-setup-and-configuration)):

```bash
seqmat setup
```

This downloads prebuilt data (genes.db + FASTA) from SeqMat's S3 bucket into a dedicated directory on your machine. No config file is required: the installer can add `SEQMAT_DATA_DIR` to your shell config so SeqMat finds the data in future sessions.

### From Source
```bash
git clone https://github.com/yourusername/seqmat.git
cd seqmat
pip install -e .
```

## Quick Start

### 1. Basic Sequence Operations

```python
from seqmat import SeqMat

# Create a sequence
seq = SeqMat("ATCGATCGATCG", name="my_sequence")
print(f"Length: {len(seq)}")  # Length: 12
print(f"Sequence: {seq.seq}")  # Sequence: ATCGATCGATCG

# Apply mutations
seq.apply_mutations([
    (3, "C", "G"),       # SNP: C->G at position 3
    (6, "-", "AAA"),     # Insertion: insert AAA at position 6  
    (10, "TC", "-")      # Deletion: delete TC at position 10
])

print(f"Mutated: {seq.seq}")
print(f"Mutations: {len(seq.mutations)}")
```

### 2. Working with Genomic Coordinates

```python
import numpy as np

# Create sequence with genomic coordinates
indices = np.arange(1000, 1012)  # Positions 1000-1011
seq = SeqMat("ATCGATCGATCG", indices=indices, name="chr1:1000-1011")

# Slice by genomic position
subseq = seq[1003:1008]  # Extract positions 1003-1007
print(f"Subsequence: {subseq.seq}")

# Access single positions
base = seq[1005]  # Get base at position 1005
print(f"Base at 1005: {base['nt'].decode()}")
```

### 3. Advanced Mutation Operations

```python
# Multiple mutations with validation
mutations = [
    (1002, "T", "A"),     # SNP
    (1005, "-", "GGG"),   # Insertion
    (1008, "GAT", "-"),   # Deletion
    (1003, "CG", "AT")    # Complex substitution
]

seq.apply_mutations(mutations)

# Mutation history
for mut in seq.mutations:
    print(f"{mut['type']} at {mut['pos']}: {mut['ref']} -> {mut['alt']}")

# Check for conflicts (automatically validated)
conflicting_muts = [
    (1003, "C", "T"),
    (1003, "G", "A")  # Overlaps with previous
]
# This will warn about conflicts and skip invalid mutations
```

### 4. Sequence Transformations

```python
# Complement and reverse complement
complement = seq.complement()
rev_comp = seq.clone()
rev_comp.reverse_complement()

print(f"Original:    {seq.seq}")
print(f"Complement:  {complement.seq}")
print(f"Rev comp:    {rev_comp.seq}")

# Remove regions (e.g., splice out introns)
introns = [(1003, 1005), (1008, 1009)]
spliced = seq.remove_regions(introns)
print(f"Spliced:     {spliced.seq}")
```

### 5. Working with Genes (requires genomics data setup)

```python
from seqmat import Gene, setup_genomics_data

# One-time setup: downloads prebuilt genes.db + FASTA from SeqMat's S3 bucket (no rebuild)
# Uses a dedicated data dir on your machine (default: platform user data dir)
setup_genomics_data("/path/to/data", organism="hg38")  # or omit path to use default location

# Load a gene
kras = Gene.from_file("KRAS", organism="hg38")
print(kras)  # Gene: KRAS, ID: ENSG00000133703, Chr: 12, Transcripts: 8

# Access transcripts
primary = kras.transcript()  # Primary transcript
all_transcripts = list(kras)  # All transcripts

# Analyze gene structure
acceptors, donors = kras.splice_sites()
print(f"Unique splice sites: {len(acceptors)} acceptors, {len(donors)} donors")
```

### 6. Transcript Analysis

```python
# Get a transcript
transcript = kras.transcript()
print(f"Transcript: {transcript.transcript_id}")
print(f"Protein coding: {transcript.protein_coding}")

# Access exon/intron structure
print(f"Exons: {len(transcript.exons)}")
print(f"Introns: {len(transcript.introns)}")

# Generate mature mRNA (spliced)
transcript.generate_mature_mrna()
print(f"Mature mRNA length: {len(transcript.mature_mrna)} bp")

# Generate protein (if protein-coding)
if transcript.protein_coding:
    protein = transcript.generate_protein()
    print(f"Protein length: {len(transcript.protein)} aa")
```

### 7. Loading from FASTA Files

```python
# Direct FASTA loading
genomic_seq = SeqMat.from_fasta_file(
    "chr12.fasta",
    "chr12", 
    start=25398284,
    end=25398384
)

# Apply mutations to genomic sequence
genomic_seq.apply_mutations([
    (25398290, "G", "A"),  # Pathogenic variant
    (25398300, "-", "T")   # Novel insertion
])
```

## Data Setup and Configuration

**Typical workflow:** `pip install seqmat` → run `seqmat setup` once (or `setup_genomics_data()` with optional path) → data is downloaded to a dedicated directory; SeqMat can add `SEQMAT_DATA_DIR` to your shell so future use is config-less. No manual config or rebuild needed.

SeqMat uses a flexible configuration system. When `SEQMAT_DATA_DIR` is set, all paths are derived from it (config-less). Otherwise paths are read from a config file (written by setup if needed).

### Quick Setup for Supported Organisms

For ordinary **pip install** users, setup **downloads prebuilt data** (genes.db and FASTA) from the SeqMat S3 bucket into a dedicated location on your machine. **No rebuild is required.** Default data location is the platform user data dir (e.g. `~/Library/Application Support/seqmat` on macOS). For the built-in organisms (hg38, mm39):

```python
from seqmat import setup_genomics_data

# Default: download prebuilt genes.db + FASTA from S3 (fast)
setup_genomics_data("/path/to/your/data", organism="hg38")
setup_genomics_data("/path/to/your/data", organism="mm39")

# Optional: build from sources (GTF + FASTA) instead of using prebuilt data
setup_genomics_data("/path/to/your/data", organism="hg38", from_prebuilt=False, n_jobs=8)
```

Or via command line:
```bash
# Default: download prebuilt data from S3
seqmat setup --path /path/to/your/data --organism hg38

# Build from GTF/FASTA sources instead
seqmat setup --path /path/to/your/data --organism hg38 --build-from-sources
```

### Configuration System

**Config-less mode (recommended):**  
When you run `seqmat setup` (or `setup_genomics_data()`), SeqMat can append to your shell config (e.g. `.zshrc` or `.bashrc`):

```bash
export SEQMAT_DATA_DIR=/path/to/seqmat_base
export SEQMAT_DEFAULT_ORGANISM=hg38
```

After you open a new terminal (or run `source ~/.zshrc`), SeqMat uses only these environment variables. All paths are derived: `{SEQMAT_DATA_DIR}/{organism}/genes.db`, `{organism}.fa`, etc. You do not need to create or edit any config file.

**Config file mode** (when `SEQMAT_DATA_DIR` is not set):

1. **`SEQMAT_CONFIG_FILE`** – exact path to config.json.
2. **`SEQMAT_CONFIG_DIR`** – directory; config is `{SEQMAT_CONFIG_DIR}/config.json`.
3. **Default** – platform config dir (e.g. `~/.config/seqmat/config.json` on Linux).

To see which file is used: `from seqmat import get_config_file; print(get_config_file())`.  
To see if config-less: `from seqmat import get_data_base; print(get_data_base())`.

The config file (when used) contains:

- **Organism paths**: Where each organism's data is stored
- **Default organism**: Which organism to use when none specified
- **Directory structure**: Customizable folder names
- **Data source URLs**: Where to download reference data

#### Configuration File Structure

```json
{
  "default_organism": "hg38",
  "directory_structure": {
    "chromosomes": "chromosomes",
    "annotations": "annotations"
  },
  "hg38": {
    "BASE": "/path/to/data/hg38",
    "CHROM_SOURCE": "/path/to/data/hg38",
    "MRNA_PATH": "/path/to/data/hg38",
    "genes_db": "/path/to/data/hg38/genes.db",
    "fasta_full_genome": "/path/to/data/hg38/hg38.fa"
  },
  "mm39": {
    "BASE": "/path/to/data/mm39",
    "CHROM_SOURCE": "/path/to/data/mm39",
    "MRNA_PATH": "/path/to/data/mm39",
    "genes_db": "/path/to/data/mm39/genes.db",
    "fasta_full_genome": "/path/to/data/mm39/mm39.fa"
  }
}
```

### Directory Structure

When you run `setup_genomics_data`, each organism directory contains a single **genes.db** (SQLite) and a full-genome FASTA:

```
/your/data/path/
├── hg38/                          # Organism directory
│   ├── genes.db                   # Gene annotations (SQLite; from prebuilt or built from GTF)
│   ├── hg38.fa                    # Full genome FASTA
│   ├── temp/                      # Temporary files (when building from sources)
│   ├── missplicing/
│   └── oncosplice/
└── mm39/                          # Mouse (same layout)
    ├── genes.db
    ├── mm39.fa
    └── ...
```

Prebuilt data is downloaded from SeqMat's S3 bucket (default: `seqmat-prebuilt-public`). To upload or update prebuilt assets, use layout: `{bucket}/{organism}/genes.db`, `{bucket}/{organism}/{organism}.fa.gz`, and optionally `{bucket}/{organism}/conservation.db` or `conservation.pkl`.

### Adding Support for New Organisms

To add a new organism (e.g., dm6 for Drosophila), you have several options:

#### Option 1: Modify Configuration (Recommended)

1. **Update the default organism data** by editing your config or using the API:

```python
from seqmat.config import load_config, save_config

# Load current config
config = load_config()

# Add new organism with data source URLs
config['dm6'] = {
    'name': 'Drosophila melanogaster (Fruit fly)',
    'urls': {
        'fasta': 'https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz',
        'gtf': 'https://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf.gz'
    }
}

# Save updated config
save_config(config)

# Now you can use it
setup_genomics_data("/path/to/data", organism="dm6")
```

#### Option 2: Extend Default Data Sources

For permanent additions, modify `seqmat/config.py`:

```python
# In seqmat/config.py, add to DEFAULT_ORGANISM_DATA:
DEFAULT_ORGANISM_DATA = {
    'hg38': { ... },
    'mm39': { ... },
    'dm6': {
        'name': 'Drosophila melanogaster (Fruit fly)',
        'urls': {
            'fasta': 'https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz',
            'gtf': 'https://ftp.ensembl.org/pub/release-109/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.32.109.gtf.gz'
        }
    }
}
```

#### Option 3: Manual Setup for Custom Data

For completely custom data: create the organism directory, place `genes.db` (SQLite, same schema as built by setup) and a full-genome FASTA (e.g. `custom.fa`), then point config at that directory:

```python
from pathlib import Path
from seqmat.config import load_config, save_config

base_path = Path("/path/to/data/custom_genome")
base_path.mkdir(parents=True, exist_ok=True)
# Place genes.db and e.g. custom_genome.fa in base_path, then:

config = load_config()
config['custom_genome'] = {
    'BASE': str(base_path),
    'CHROM_SOURCE': str(base_path),
    'MRNA_PATH': str(base_path),
    'genes_db': str(base_path / 'genes.db'),
    'fasta_full_genome': str(base_path / 'custom_genome.fa')
}
save_config(config)
```

### Data Sources

**Default (prebuilt download):**  
SeqMat downloads prebuilt `genes.db` and FASTA from the SeqMat S3 bucket (default: `seqmat-prebuilt-public`). No rebuild; data is ready to use.

**Build-from-sources** (when you pass `from_prebuilt=False` or `--build-from-sources`):
- GTF: Ensembl (e.g. release-111 for hg38, release-112 for mm39)
- FASTA: UCSC Genome Browser (hg38, mm39)
- Conservation: same prebuilt bucket, per organism (`conservation.db` or `.pkl`)
- Expression (human only): GTEx v8 median TPM data

### Overriding the prebuilt data source

By default, prebuilt data is downloaded from SeqMat's S3 bucket. To use a different bucket or mirror, set one of:

- **Environment:** `export SEQMAT_PREBUILT_DATA_BASE_URL=https://your-bucket.s3.region.amazonaws.com`
- **Config file:** add `"prebuilt_data_base_url": "https://..."` to your seqmat config (see `get_config_file()`).

### Configuration Management

**View current configuration:**
```python
from seqmat.config import load_config, get_available_organisms, get_default_organism

print("Default organism:", get_default_organism())
print("Available organisms:", get_available_organisms())
print("Full config:", load_config())
```

**Reset to defaults:**
```python
from seqmat.config import DEFAULT_SETTINGS, save_config
save_config(DEFAULT_SETTINGS.copy())
```

**Change default organism:**
```python
from seqmat.config import load_config, save_config

config = load_config()
config['default_organism'] = 'mm39'  # Switch default to mouse
save_config(config)
```

### Rebuilding data from sources

If you need to regenerate `genes.db` from a GTF (e.g. a different Ensembl release or custom annotations), use build-from-sources. This downloads the GTF and FASTA from Ensembl/UCSC, fetches conservation from the prebuilt bucket, parses the GTF, and writes `genes.db` locally.

- **CLI:** `seqmat setup --path /path/to/data --organism hg38 --build-from-sources`
- **Python:** `setup_genomics_data("/path/to/data", organism="hg38", from_prebuilt=False, n_jobs=8)`

### Troubleshooting Setup

**Common issues:**

1. **"Organism not configured"** - Run setup first or check config
2. **Download failures** - Check internet connection and URLs
3. **Permission errors** - Ensure write access to data directory
4. **Disk space** - Human genome data requires ~4GB, mouse ~3GB

**Debugging:**
```python
# Check what's configured
from seqmat import list_available_organisms, print_data_summary

print("Configured organisms:", list_available_organisms()) 
print_data_summary()  # Detailed status

# Verify file paths
from seqmat.config import get_organism_config
config = get_organism_config("hg38")
print("Paths:", config)
```

### Configuration Best Practices

**For single-user systems:**
- Use the default location (see resolution order above) or set `SEQMAT_CONFIG_FILE` or `SEQMAT_CONFIG_DIR`.
- Set up organisms as needed with `setup_genomics_data()`

**For multi-user or shared systems:**
- Create a shared data directory: `/shared/genomics_data/`
- Point all users' configs to the same data paths
- Consider using environment variables for paths

**For ephemeral environments (e.g. Run.ai, Docker):**
- Store `config.json` in persistent storage and point to it: `export SEQMAT_CONFIG_FILE=/path/to/config.json` or `export SEQMAT_CONFIG_DIR=/path/to/dir` (config will be `dir/config.json`).

**For development/testing:**
- Use separate config files or directories
- Set `default_organism` to your most-used organism
- Keep test data in separate locations

**Configuration sharing:**
```python
# Export configuration to share with others
from seqmat.config import load_config
import json

config = load_config()
with open('seqmat_config_template.json', 'w') as f:
    json.dump(config, f, indent=2)
```

When building from sources, setup produces per organism:
- **genes.db** (SQLite) in the organism directory
- **Full-genome FASTA** (e.g. `hg38.fa`) in the organism directory
- **Expression data** (human only) merged into gene records

## Data Inspection and Management

### Python API

Once data is set up, you can inspect what's available:

```python
from seqmat import (
    list_supported_organisms, list_available_organisms,
    get_organism_info, list_gene_biotypes, count_genes,
    get_gene_list, search_genes, print_data_summary
)

# Check supported and configured organisms
print("Supported:", list_supported_organisms())  # ['hg38', 'mm39']  
print("Configured:", list_available_organisms())  # Organisms with data

# Get detailed organism information
info = get_organism_info('hg38')
print(f"Gene types: {info['data_available']['biotypes']}")
print(f"Chromosomes: {len(info['data_available']['chromosomes'])}")

# Explore gene biotypes and counts
biotypes = list_gene_biotypes('hg38')  # ['protein_coding', 'lncRNA', ...]
counts = count_genes('hg38')  # {'protein_coding': 19234, 'lncRNA': 7805, ...}

# Get gene lists
protein_genes = get_gene_list('hg38', 'protein_coding', limit=10)
print(f"First 10 protein-coding genes: {protein_genes}")

# Search for specific genes
results = search_genes('hg38', 'KRAS')  # Find genes matching 'KRAS'
kras_genes = search_genes('hg38', 'K', biotype='protein_coding', limit=5)

# Print comprehensive summary
print_data_summary()  # Formatted overview of all data
```

### Command Line Interface

SeqMat provides a comprehensive CLI for data management. The CLI automatically detects available organisms from your configuration:

```bash
# Default: download prebuilt data (no rebuild)
seqmat setup
seqmat setup --path /your/data/path --organism hg38

# Build from GTF/FASTA instead (for custom or updated annotations)
seqmat setup --path /your/data/path --organism hg38 --build-from-sources

seqmat setup --help  # Shows path default and --organism choices

# Check what organisms are supported/configured
seqmat organisms

# Get comprehensive data summary
seqmat summary

# List gene biotypes for an organism  
seqmat biotypes --organism hg38

# Count genes by biotype
seqmat count --organism hg38                    # All biotypes
seqmat count --organism hg38 --biotype protein_coding  # Specific biotype

# List genes
seqmat list --organism hg38 --biotype protein_coding --limit 20

# Search for genes by name
seqmat search --organism hg38 --query KRAS
seqmat search --organism hg38 --query K --biotype protein_coding --limit 10

# Get detailed organism information
seqmat info --organism hg38
```

### Example CLI Output

```bash
$ seqmat summary
🧬 SeqMat Genomics Data Summary
========================================
📊 Total: 2 organisms, 15 biotypes, 47,832 genes

🌍 Supported Organisms:
   hg38: Homo sapiens (Human) - ✅ Configured
   mm39: Mus musculus (Mouse) - ✅ Configured

📁 HG38 Data:
   Gene Types:
     protein_coding: 19,234 genes
     lncRNA: 7,805 genes
     pseudogene: 14,723 genes
     ...
   Chromosomes: 25 available (chr1, chr2, chr3, chr4, chr5...)

📁 MM39 Data:
   Gene Types:
     protein_coding: 21,815 genes
     lncRNA: 8,032 genes
     ...
```

## API Reference

### SeqMat Class

**Core Methods:**
- `SeqMat(nucleotides, indices=None, name="wild_type")`: Create sequence
- `SeqMat.from_fasta_file(path, chrom, start, end)`: Load from FASTA
- `apply_mutations(mutations)`: Apply SNPs/indels
- `clone(start=None, end=None)`: Create copy
- `remove_regions(regions)`: Remove specified intervals
- `complement()`: Get complement sequence
- `reverse_complement()`: Reverse complement in place

**Properties:**
- `seq`: Current sequence string
- `reference_seq`: Original reference sequence
- `index`: Genomic coordinates
- `mutations`: List of applied mutations
- `mutated_positions`: Set of mutated positions

### Gene Class

- `Gene.from_file(gene_name, organism=None)`: Load gene from database (uses default organism if None)
- `transcript(tid=None)`: Get transcript by ID or primary
- `splice_sites()`: Get all splice site positions
- `primary_transcript`: Primary transcript ID

### Transcript Class

- `generate_pre_mrna()`: Create pre-mRNA sequence
- `generate_mature_mrna()`: Create spliced mRNA
- `generate_protein()`: Translate to protein
- `exons`/`introns`: Genomic coordinates
- `protein_coding`: Boolean flag

### Configuration Management

**Core Functions:**
- `load_config()`: Load configuration (file path from `get_config_file()`)
- `get_config_file()`: Return the active config file path
- `get_data_base()`: Return `SEQMAT_DATA_DIR` path if set (config-less mode), else None
- `save_config(config)`: Save configuration to file
- `get_default_organism()`: Get default organism (from env in config-less mode, else config)
- `get_available_organisms()`: Get list of configured organisms (from data dir in config-less mode)
- `get_organism_config(organism=None)`: Get file paths for organism (uses default if None)
- `get_organism_info(organism)`: Get organism metadata
- `get_prebuilt_data_base_url()`: Base URL for prebuilt downloads (env > config > default bucket)

**Configuration Examples:**
```python
from seqmat.config import *

# Check current settings
print("Default:", get_default_organism())         # 'hg38'
print("Available:", get_available_organisms())    # ['hg38', 'mm39']  
config = get_organism_config('hg38')              # File paths

# Modify settings
config = load_config()
config['default_organism'] = 'mm39'
save_config(config)
```

### Data Inspection Utilities

**Organism Management:**
- `list_supported_organisms()`: Get all supported organisms (dynamic from config)
- `list_available_organisms()`: Get configured organisms  
- `get_organism_info(organism)`: Detailed organism information
- `setup_genomics_data(basepath, organism=None, force=False, pickup=False, n_jobs=None, from_prebuilt=True)`: Set up data. By default downloads prebuilt genes.db and FASTA from S3; use `from_prebuilt=False` to build from GTF/FASTA. `n_jobs` controls GTF parsing workers when building from sources.

**Gene Discovery:**
- `list_gene_biotypes(organism)`: Get available gene types
- `count_genes(organism, biotype=None)`: Count genes by type
- `get_gene_list(organism, biotype, limit=None)`: List gene names  
- `search_genes(organism, query, biotype=None, limit=10)`: Search genes by name

**Data Summary:**
- `data_summary()`: Complete data overview (programmatic)
- `print_data_summary()`: Formatted data summary (human-readable)

### Command Line Interface

**Setup Commands:**
- `seqmat setup --path PATH --organism {dynamic_list} [--force] [--build-from-sources]`: Default: download prebuilt data from S3; use `--build-from-sources` to build genes.db from GTF/FASTA
- `seqmat organisms`: List organism status and availability

**Exploration Commands:**  
- `seqmat summary`: Data overview for all configured organisms
- `seqmat info --organism ORG`: Detailed organism info
- `seqmat biotypes --organism ORG`: List gene biotypes
- `seqmat count --organism ORG [--biotype TYPE]`: Count genes
- `seqmat list --organism ORG --biotype TYPE [--limit N]`: List genes
- `seqmat search --organism ORG --query PATTERN [--biotype TYPE] [--limit N]`: Search genes

**Note**: All CLI commands automatically use your configured default organism when `--organism` is omitted, and available organism choices are dynamically loaded from your configuration.

## Performance

SeqMat is optimized for performance:

- **Vectorized operations**: NumPy-based sequence operations
- **Memory efficient**: Structured arrays for sequence storage
- **Fast slicing**: O(1) genomic coordinate access
- **Conflict detection**: Efficient mutation validation
- **Lazy loading**: Sequences loaded on demand

## Dependencies

**Core (always required):**
- numpy >= 1.20.0, <3
- pandas >= 2.2.0 (required for NumPy 2 compatibility)
- pysam >= 0.19.0
- requests >= 2.26.0
- tqdm >= 4.62.0

**Optional:**
- gtfparse >= 1.2.0 (for genomics data setup)
- pyarrow >= 14 (included; used for parquet export in `save_seqmat`)

## Examples

Check the `examples/` directory for:
- Basic sequence manipulation
- Mutation analysis workflows
- Gene structure analysis
- Comparative genomics
- Performance benchmarks

## Contributing

Contributions welcome! Please see CONTRIBUTING.md for guidelines.

## License

MIT License - see LICENSE file for details.

## Citation

If you use SeqMat in your research:

```
SeqMat: Lightning-fast genomic sequence matrix library
[Your Name], 2024
GitHub: https://github.com/yourusername/seqmat
```

## Support

- **Documentation**: https://seqmat.readthedocs.io
- **Issues**: https://github.com/yourusername/seqmat/issues
- **Discussions**: https://github.com/yourusername/seqmat/discussions