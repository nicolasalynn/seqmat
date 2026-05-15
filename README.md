<p align="center">
  <img src="docs/assets/logo.svg" alt="SeqMat" width="540">
</p>

[![PyPI](https://img.shields.io/pypi/v/seqmat.svg)](https://pypi.org/project/seqmat/)
[![Python](https://img.shields.io/pypi/pyversions/seqmat.svg)](https://pypi.org/project/seqmat/)
[![Tests](https://github.com/nicolasalynn/seqmat/actions/workflows/test.yml/badge.svg)](https://github.com/nicolasalynn/seqmat/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/nicolasalynn/seqmat/branch/main/graph/badge.svg)](https://codecov.io/gh/nicolasalynn/seqmat)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Downloads](https://static.pepy.tech/badge/seqmat/month)](https://pepy.tech/project/seqmat)

A Python library for genomic sequences that keeps genomic coordinates and mutation history attached to the bases.

SeqMat stores a DNA sequence as a NumPy structured array — `(nt, ref, index, mut_type, valid)` in parallel columns rather than a byte string. Slicing, mutation, complement, and splicing preserve genomic coordinates and the reference state through every transform, and the mutation history is recorded as a side product. The library also includes a gene/transcript model loaded from a single SQLite file built from Ensembl annotations, and a position → gene lookup over a sidecar index.

```python
from seqmat import Gene, SeqMat

# What gene is at chr12:25,245,350?
Gene.from_position("12", 25_245_350)            # [Gene(KRAS)]

# Load it, assemble the mature mRNA, translate.
kras = Gene.from_file("KRAS")
tx   = kras.transcript()
tx.generate_mature_mrna()
tx.generate_protein()

# Apply mutations; history is kept on the SeqMat.
seq = SeqMat("ATCGATCGATCG")
seq.apply_mutations([(3, "C", "G"), (6, "-", "AAA"), (10, "TC", "-")])
seq.mutations
```

Worked examples:

- [KRAS G12D — coordinate to protein](examples/kras_g12d_analysis.ipynb): chromosome coordinate → mature mRNA → protein, comparing wild type to the G12D mutant.
- [SeqMat vs Biopython, side by side](examples/seqmat_vs_biopython.ipynb): the same tasks in both libraries, with timings, calling out where each one wins.

---

## Install

```bash
pip install seqmat
seqmat setup                                    # one-time: downloads hg38 genes.db + FASTA (~4 GB)
```

`seqmat setup` appends `SEQMAT_DATA_DIR` and `SEQMAT_DEFAULT_ORGANISM` to your shell rc so the data is found automatically next session. Mouse: `seqmat setup --organism mm39`.

If you only want `SeqMat` for in-memory sequence work, you can skip `seqmat setup` — `Gene` and `Transcript` are the only objects that need the gene database.

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

### Position → gene (added in 1.4.0)

```python
from seqmat import Gene, gene_names_at_position

Gene.from_position("12", 25_245_350)            # [Gene(KRAS)] — point query
Gene.from_position("chr12", (25_200_000, 25_300_000))  # all overlapping genes in a range
gene_names_at_position("X", 100_000)            # names only (no BLOB load)
```

Backed by a per-chromosome sorted NumPy index persisted as a sidecar `gene_locations.npz` next to `genes.db`. Built lazily on first call; fresh `seqmat setup` builds also emit it.

### Loading from FASTA

```python
seq = SeqMat.from_fasta_file("chr12.fasta", "chr12", start=25398284, end=25398384)
seq.apply_mutations([(25398290, "G", "A")])     # G12D, a common KRAS variant
```

## Performance

Measured on an M-series Mac, hg38, one core, warm caches. Reproducible via the scripts in `benchmarks/`.

### Headline numbers

| Operation                                          | Time     |
| -------------------------------------------------- | -------: |
| `gene_names_at_position(chrm, pos)`                | 2.5 µs   |
| KRAS mature mRNA assembly (splice + translate)     | 0.2 ms   |
| 1,000-SNP batch on 45 kb sequence (with history)   | 0.5 ms   |
| `Gene.from_file("KRAS")` (SQLite + unpickle)       | 24 ms    |
| `Gene.from_position(chrm, pos)` end-to-end         | 24 ms    |

### Position → gene

Same 63,241 hg38 gene intervals, same 10,000 random point queries. Different libraries are designed for different access patterns, so we report two workloads. Reproduce with `python benchmarks/bench_position_lookup.py`.

Per-query (one coordinate, one answer, in a loop — the `Gene.from_position` pattern):

| Implementation                          | Per query |
| --------------------------------------- | --------: |
| SeqMat (locator)                        | 2.5 µs    |
| Python `dict` + `bisect`                | 21 µs     |
| pandas (`groupby` chrm + boolean mask)  | 79 µs     |
| PyRanges constructed per call           | 2.0 ms    |

Batched (the whole query set in one call — PyRanges' native idiom):

| Implementation                | Per query |
| ----------------------------- | --------: |
| PyRanges `.join`              | 2.07 µs   |
| SeqMat (serial loop)          | 2.59 µs   |

PyRanges and SeqMat are in roughly the same range when each is used the way it was designed for. Constructing a `PyRanges` per call is the slow path and not the PyRanges authors' intent — included here only because it's a common mistake.

### Sequence operations vs Biopython

50 kb sequence, 1,000 SNPs / 10 introns; 1 Mb reverse-complement. Reproduce with `python benchmarks/bench_sequence_ops.py`.

| Workload                                       | Biopython / native              | SeqMat  |
| ---------------------------------------------- | ------------------------------: | ------: |
| 1,000 SNPs (no history)                        | `MutableSeq`: 0.21 ms           | —       |
| 1,000 SNPs with full mutation history          | —                               | 0.5 ms  |
| 10-intron splice (50 kb)                       | `str.join`: 0.002 ms            | 0.8 ms  |
| Reverse-complement (1 Mb)                      | `Seq.reverse_complement`: 0.57 ms | 9.4 ms |

For raw byte-string throughput, `Bio.Seq` and `bytearray` are faster. SeqMat's structured array adds bookkeeping that those types don't carry — genomic coordinates that survive indels, the reference held alongside the current sequence, and a mutation history list. The trade-off makes sense when you need that bookkeeping anyway (the KRAS G12D walkthrough is a typical case) and doesn't when you don't.

### Implementation notes

NumPy structured arrays, LUT-based complement, `np.searchsorted` on sorted exon and gene-start coordinates, FASTA range-scoped reads, copy-on-write `clone()`. See [`seqmat/seqmat.py`](seqmat/seqmat.py) and [`seqmat/locator.py`](seqmat/locator.py).

## Command-line interface

```bash
seqmat setup [--organism hg38|mm39] [--path PATH] [--build-from-sources]
seqmat summary
seqmat info --organism hg38
seqmat search --organism hg38 --query KRAS
seqmat list   --organism hg38 --biotype protein_coding --limit 20
seqmat count  --organism hg38
```

## Data setup

By default `seqmat setup` downloads prebuilt `genes.db` and FASTA from the SeqMat S3 bucket. To regenerate `genes.db` from a specific Ensembl release or a custom GTF:

```bash
seqmat setup --organism hg38 --build-from-sources
```

Custom organisms, mirror buckets, ephemeral environments (Docker / Run.ai), shared multi-user installs, and the full configuration system are documented in [docs/SETUP.md](docs/SETUP.md).

## API at a glance

```python
from seqmat import (
    SeqMat,                           # coordinate-tracked sequence with mutation history
    Gene, Transcript,                 # gene/transcript model
    gene_names_at_position,           # name-only positional lookup
    build_location_index,             # force-rebuild the position index
    setup_genomics_data,              # programmatic setup
    search_genes, available_genes,    # discovery helpers
)
```

Key classes:

- `SeqMat` — `apply_mutations`, `clone`, `complement`, `reverse_complement`, `remove_regions`, `from_fasta_file`
- `Gene` — `from_file`, `from_position`, `get`, `transcript`, `splice_sites`, `primary_transcript`
- `Transcript` — `generate_pre_mrna`, `generate_mature_mrna`, `generate_protein`, `exons`, `introns`

Full API reference (auto-generated from docstrings): [nicolasalynn.github.io/seqmat](https://nicolasalynn.github.io/seqmat/seqmat.html).

## Requirements

Python ≥ 3.10. Core deps: `numpy`, `pandas`, `pyarrow`, `pysam`, `requests`, `tqdm`, `platformdirs`. Optional: `lmdb` for faster gene loading on large workloads (`pip install seqmat[lmdb]`).

## A note on indexing

`SeqMat` defaults to 1-based indices to match genomic conventions (UCSC, GenBank). If you'd rather work in 0-based offsets — e.g. for a comparison with `bytearray` or `Bio.Seq` — pass `indices=np.arange(len(seq))` at construction.

## Contributing

PRs welcome. Run the test suite with `pytest tests/`. Comparative benchmarks are under `benchmarks/`.

## License

MIT — see [LICENSE](LICENSE).

## Citation

If SeqMat is useful in your research:

```
Lynn Vila, N. (2025). SeqMat: a genomic sequence library with
coordinate tracking and mutation history.
https://github.com/nicolasalynn/seqmat
```
