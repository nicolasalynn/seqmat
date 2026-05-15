# SeqMat

[![PyPI](https://img.shields.io/pypi/v/seqmat.svg)](https://pypi.org/project/seqmat/)
[![Python](https://img.shields.io/pypi/pyversions/seqmat.svg)](https://pypi.org/project/seqmat/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Downloads](https://static.pepy.tech/badge/seqmat/month)](https://pepy.tech/project/seqmat)

**Fast, vectorized genomic sequences with first-class mutation tracking.**

SeqMat treats a DNA sequence as a NumPy-backed matrix of (nucleotide, genomic position) — so slicing, mutation, complement, and splicing are all vectorized array operations. It ships with a compact gene/transcript model that loads from a single indexed SQLite file built from Ensembl annotations, plus a position → gene lookup that resolves overlapping genes in microseconds.

```python
from seqmat import Gene, SeqMat

# What gene is at chr12:25,245,350?
Gene.from_position("12", 25_245_350)            # [Gene(KRAS)]

# Load it, assemble the mature mRNA, translate.
kras = Gene.from_file("KRAS")
tx   = kras.transcript()
tx.generate_mature_mrna()                       # 0.2 ms
tx.generate_protein()

# Mutate a sequence with full history and conflict detection.
seq = SeqMat("ATCGATCGATCG")
seq.apply_mutations([(3, "C", "G"), (6, "-", "AAA"), (10, "TC", "-")])
seq.mutations                                   # [(SNP, 3, C, G), (INS, 6, -, AAA), (DEL, 10, TC, -)]
```

---

## Install

```bash
pip install seqmat
seqmat setup                                    # one-time: downloads hg38 genes.db + FASTA (~4 GB)
```

`seqmat setup` writes a small env block (`SEQMAT_DATA_DIR`, `SEQMAT_DEFAULT_ORGANISM`) to your shell rc so the data is found automatically next session. Mouse: `seqmat setup --organism mm39`. No config file needed.

> If you only want `SeqMat` for in-memory sequence work, you can skip `seqmat setup` entirely — `Gene` / `Transcript` are the only things that need the gene database.

## Quick start

### Sequence operations

```python
from seqmat import SeqMat
import numpy as np

seq = SeqMat("ATCGATCGATCG", indices=np.arange(1000, 1012))
seq[1005]                                       # base at genomic position 1005
seq[1003:1008].seq                              # "GATCG"
seq.reverse_complement()                        # in place
seq.remove_regions([(1003, 1005), (1008, 1009)])  # splice out introns
```

### Genes and transcripts

```python
from seqmat import Gene

kras = Gene.from_file("KRAS")
kras                                            # Gene: KRAS, ID: ENSG00000133703, Chr: 12, Transcripts: 14

tx = kras.transcript()                          # primary transcript
tx.generate_mature_mrna()
tx.generate_protein()
tx.protein[:20]                                 # 'MTEYKLVVVGAGGVGKSALT'

acceptors, donors = kras.splice_sites()         # Counter across all transcripts
```

### Position → gene (new in 1.4.0)

```python
from seqmat import Gene, gene_names_at_position

Gene.from_position("12", 25_245_350)            # [Gene(KRAS)] — point query
Gene.from_position("chr12", (25_200_000, 25_300_000))  # all overlapping genes in a range
gene_names_at_position("X", 100_000)            # names only (no BLOB load) — ~17 us
```

Backed by a per-chromosome sorted NumPy index persisted as a sidecar `gene_locations.npz` next to `genes.db`. Built lazily on first call; fresh `seqmat setup` builds also emit it for free.

### Loading from FASTA

```python
seq = SeqMat.from_fasta_file("chr12.fasta", "chr12", start=25398284, end=25398384)
seq.apply_mutations([(25398290, "G", "A")])     # G12D, the most famous KRAS variant
```

## Performance

Numbers from an M-series Mac on hg38 (one core, warm caches):

| Operation                              | Time     |
| -------------------------------------- | -------: |
| `gene_names_at_position(chrm, pos)`    | **17 µs** |
| `Gene.from_file("KRAS")` (SQLite + unpickle) | 24 ms |
| `Gene.from_position(chrm, pos)` end-to-end   | 24 ms |
| KRAS mature mRNA assembly              | 0.2 ms |
| 1,000-SNP batch on 4 kb sequence       | 19 ms |

The hot paths use NumPy structured arrays, LUT-based complement, `np.searchsorted` on sorted starts, FASTA range-scoped reads, and copy-on-write `clone()`. See `seqmat/seqmat.py` and `seqmat/locator.py`.

## Command-line interface

```bash
seqmat setup [--organism hg38|mm39] [--path PATH] [--build-from-sources]
seqmat summary                                  # what's installed
seqmat info --organism hg38
seqmat search --organism hg38 --query KRAS
seqmat list   --organism hg38 --biotype protein_coding --limit 20
seqmat count  --organism hg38
```

## Data setup

By default `seqmat setup` downloads prebuilt `genes.db` and FASTA from the SeqMat S3 bucket — no build step. To regenerate `genes.db` from a specific Ensembl release or custom GTF:

```bash
seqmat setup --organism hg38 --build-from-sources
```

For custom organisms, mirroring the prebuilt bucket, ephemeral environments (Docker / Run.ai), shared multi-user installs, and the full configuration system — see **[docs/SETUP.md](docs/SETUP.md)**.

## API at a glance

```python
from seqmat import (
    SeqMat,                           # vectorized sequence with mutation tracking
    Gene, Transcript,                 # gene/transcript model
    gene_names_at_position,           # fast name-only positional lookup
    build_location_index,             # force-rebuild the position index
    setup_genomics_data,              # programmatic setup
    search_genes, available_genes,    # discovery helpers
)
```

Key classes:

- **`SeqMat`** — `apply_mutations`, `clone`, `complement`, `reverse_complement`, `remove_regions`, `from_fasta_file`
- **`Gene`** — `from_file`, `from_position`, `transcript`, `splice_sites`, `primary_transcript`
- **`Transcript`** — `generate_pre_mrna`, `generate_mature_mrna`, `generate_protein`, `exons`, `introns`

## Requirements

Python ≥ 3.10. Core deps: `numpy`, `pandas`, `pyarrow`, `pysam`, `requests`, `tqdm`, `platformdirs`. Optional: `lmdb` (faster gene loading on large workloads — install with `pip install seqmat[lmdb]`).

## Contributing

PRs welcome. Run the test suite with `pytest tests/`. Benchmarks live under `tests/bench_*.py`.

## License

MIT — see [LICENSE](LICENSE).

## Citation

If SeqMat is useful in your research, please cite:

```
Lynn Vila, N. (2025). SeqMat: a fast, vectorized genomic sequence library
with mutation tracking. https://github.com/nicolasalynn/seqmat
```
