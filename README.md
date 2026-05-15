<p align="center">
  <img src="docs/assets/logo.svg" alt="SeqMat" width="540">
</p>

[![PyPI](https://img.shields.io/pypi/v/seqmat.svg)](https://pypi.org/project/seqmat/)
[![Python](https://img.shields.io/pypi/pyversions/seqmat.svg)](https://pypi.org/project/seqmat/)
[![Tests](https://github.com/nicolasalynn/seqmat/actions/workflows/test.yml/badge.svg)](https://github.com/nicolasalynn/seqmat/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/nicolasalynn/seqmat/branch/main/graph/badge.svg)](https://codecov.io/gh/nicolasalynn/seqmat)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Downloads](https://static.pepy.tech/badge/seqmat/month)](https://pepy.tech/project/seqmat)

**Vectorized genomic sequences as a coordinate-tracked matrix, with first-class mutation history.**

SeqMat stores a DNA sequence as a NumPy structured array of `(nt, ref, index, mut_type, valid)` — five parallel columns rather than a byte string. That extra book-keeping is the whole point: slicing, mutation, complement, and splicing all preserve **genomic coordinates** and **mutation provenance** through every transform. You get the convenience of "what was originally at chr12:25,245,350 vs what's there now" for free; in `Bio.Seq` or `bytearray` you'd be maintaining offset tables and a side-car history dict yourself.

It ships with a compact gene/transcript model that loads from a single indexed SQLite file built from Ensembl annotations, plus a position → gene lookup that resolves overlapping genes in microseconds.

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

📘 **Walkthrough:** [Take a genomic coordinate through to the protein-level diff of a KRAS G12D oncogenic mutation](examples/kras_g12d_analysis.ipynb) — five cells, all the killer features.

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
gene_names_at_position("X", 100_000)            # names only (no BLOB load) — ~2.5 us
```

Backed by a per-chromosome sorted NumPy index persisted as a sidecar `gene_locations.npz` next to `genes.db`. Built lazily on first call; fresh `seqmat setup` builds also emit it for free.

### Loading from FASTA

```python
seq = SeqMat.from_fasta_file("chr12.fasta", "chr12", start=25398284, end=25398384)
seq.apply_mutations([(25398290, "G", "A")])     # G12D, the most famous KRAS variant
```

## Performance

### Headline numbers

M-series Mac on hg38, one core, warm caches:

| Operation                                    | Time       |
| -------------------------------------------- | ---------: |
| `gene_names_at_position(chrm, pos)`          | **2.5 µs** |
| KRAS mature mRNA assembly (splice + translate) | 0.2 ms   |
| 1,000-SNP batch on 45 kb sequence (with history) | 0.5 ms |
| `Gene.from_file("KRAS")` (SQLite + unpickle) | 24 ms      |
| `Gene.from_position(chrm, pos)` end-to-end   | 24 ms      |

### Position → gene, head-to-head

Same 63,241 hg38 gene intervals, same 10,000 random point queries, same machine. **Two workloads** — different libraries are optimized for different patterns, so we report both. Reproduce with `python benchmarks/bench_position_lookup.py`.

**Per-query (the `Gene.from_position` pattern: one coordinate, one answer, in a loop):**

| Implementation                    | Per query  | Relative      |
| --------------------------------- | ---------: | ------------: |
| **SeqMat (locator)**              | **2.5 µs** | **1×**        |
| Python `dict` + `bisect`          | 21 µs      | 8× slower    |
| pandas (`groupby` chrm + mask)    | 79 µs      | 31× slower   |
| PyRanges constructed per call     | 2.0 ms     | 800× slower *(anti-pattern)* |

**Batched (hand the whole query set to the library at once, PyRanges' native idiom):**

| Implementation                    | Per query  | Relative      |
| --------------------------------- | ---------: | ------------: |
| PyRanges `.join` (one call)       | 2.07 µs    | 1×            |
| SeqMat (serial loop)              | 2.59 µs    | 1.25× slower |

> **The honest story:** SeqMat is built for the per-query pattern and wins it convincingly — 8× over a careful hand-rolled index, 31× over pandas. PyRanges is built for batch interval-set algebra and wins that workload, but only by ~25%. Don't construct a PyRanges per call (the "PyRanges anti-pattern" row); use SeqMat for single-coordinate annotation loops or batch your queries into PyRanges. Either is fine; doing neither is what hurts.

### Sequence operations vs Biopython

50 kb sequence, 1,000 SNPs / 10 introns / 1 Mb reverse-complement. Reproduce with `python benchmarks/bench_sequence_ops.py`.

| Workload                                  | Biopython / native | SeqMat | Delta |
| ----------------------------------------- | -----------------: | -----: | ----: |
| 1,000 SNPs (no history)                   | 0.21 ms `MutableSeq` | — | — |
| **1,000 SNPs (with full mutation history)** | — | **0.5 ms** | **2.5× of `MutableSeq`, while tracking every change** |
| 10-intron splice (50 kb)                  | 0.002 ms `str.join` | 0.8 ms | 500× — Biopython has no coordinate-preserving splice |
| Reverse-complement (1 Mb)                 | 0.57 ms `Seq.reverse_complement` | 9.4 ms | 16× — reverses the 27-byte structured record array |

> **Pick the right tool:** If you just want raw byte-string throughput, `Bio.Seq` and `bytearray` will always beat a coordinate-tracked structured array. SeqMat earns the overhead when you need the *coordinates*, the *mutation history*, and the *reference vs current* state together. Workloads like the [KRAS G12D notebook](examples/kras_g12d_analysis.ipynb) — apply a few mutations at genomic coordinates, re-splice, re-translate, compare to WT — never see that overhead, because the bookkeeping is exactly what you'd have to write by hand otherwise.

### Why it's fast

NumPy structured arrays, LUT-based complement, `np.searchsorted` on sorted exon and gene-start coordinates, FASTA range-scoped reads, and copy-on-write `clone()`. See [`seqmat/seqmat.py`](seqmat/seqmat.py) and [`seqmat/locator.py`](seqmat/locator.py).

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
- **`Gene`** — `from_file`, `from_position`, `get`, `transcript`, `splice_sites`, `primary_transcript`
- **`Transcript`** — `generate_pre_mrna`, `generate_mature_mrna`, `generate_protein`, `exons`, `introns`

Full API reference: **[nicolasalynn.github.io/seqmat](https://nicolasalynn.github.io/seqmat/seqmat.html)** *(auto-generated from docstrings).*

## Requirements

Python ≥ 3.10. Core deps: `numpy`, `pandas`, `pyarrow`, `pysam`, `requests`, `tqdm`, `platformdirs`. Optional: `lmdb` (faster gene loading on large workloads — install with `pip install seqmat[lmdb]`).

## Contributing

PRs welcome. Run the test suite with `pytest tests/`. Comparative benchmarks live under `benchmarks/`.

## License

MIT — see [LICENSE](LICENSE).

## Citation

If SeqMat is useful in your research, please cite:

```
Lynn Vila, N. (2025). SeqMat: a fast, vectorized genomic sequence library
with mutation tracking. https://github.com/nicolasalynn/seqmat
```
