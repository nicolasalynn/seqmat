<p align="center">
  <img src="docs/assets/logo.svg" alt="SeqMat" width="540">
</p>

[![PyPI](https://img.shields.io/pypi/v/seqmat.svg)](https://pypi.org/project/seqmat/)
[![Python](https://img.shields.io/pypi/pyversions/seqmat.svg)](https://pypi.org/project/seqmat/)
[![Tests](https://github.com/nicolasalynn/seqmat/actions/workflows/test.yml/badge.svg)](https://github.com/nicolasalynn/seqmat/actions/workflows/test.yml)
[![codecov](https://codecov.io/gh/nicolasalynn/seqmat/branch/main/graph/badge.svg)](https://codecov.io/gh/nicolasalynn/seqmat)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Downloads](https://static.pepy.tech/badge/seqmat/month)](https://pepy.tech/project/seqmat)

A small Python library for **variant-effect and splicing workflows** — the kind where keeping a genomic coordinate and the original reference attached to each base is most of the value.

SeqMat is not a replacement for Biopython. Biopython is the right default for sequence I/O, alignments, raw byte operations, and most of bioinformatics. SeqMat fills a narrower slot: when you need to apply a variant at a genomic coordinate, splice introns out, re-translate, and ask *"what was here originally?"* — that bookkeeping lives in the data model instead of in your code.

A DNA sequence is stored as a NumPy structured array — `(nt, ref, index, mut_type, valid)` in parallel columns rather than a byte string. Slicing, mutation, complement, and splicing preserve the `index` (genomic coordinates) and the `ref` (original bases) through every transform, and mutation history is recorded as a side product. The library also includes a gene/transcript model loaded from a single SQLite file built from Ensembl annotations, and a position → gene lookup over a sidecar NumPy index.

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

## When SeqMat fits — and when it doesn't

Reach for SeqMat when:

- You have a genomic coordinate and want the gene, transcript, mature mRNA, or protein-level effect of a variant at it.
- You need to splice introns out of a pre-mRNA while keeping a coordinate map for everything that survives.
- You want a mutation history attached to the sequence without writing your own side-car dict.
- You're doing splicing or variant-effect analyses where "what was originally here?" comes up repeatedly.

Reach for Biopython (or another tool) when:

- You're reading or writing FASTA, GenBank, EMBL, alignments, BLAST output, NCBI fetch, etc. Biopython has decades of format coverage that SeqMat will not match.
- You need maximum throughput on raw byte operations — reverse-complement on long chromosomes, mass translations, big alignment passes. SeqMat's structured-array overhead is real and `Bio.Seq` will win.
- You're doing batch interval-set algebra over many feature tables — that's PyRanges' job.

The two libraries compose. Many real pipelines load FASTA via Biopython, do byte-level work in `Bio.Seq`, and reach for SeqMat where coordinate tracking and mutation history would otherwise be hand-written. See [examples/seqmat_vs_biopython.ipynb](examples/seqmat_vs_biopython.ipynb) for a head-to-head.

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

Two short statements first, then the numbers.

1. **For the operations SeqMat is built around** (coordinate-aware gene lookup, splicing-with-coordinate-preservation, mutation-tracked sequence assembly), it's fast enough that performance won't be the bottleneck in your pipeline.
2. **For raw byte operations** (reverse-complement, mass translation, simple string mutations without coordinate tracking), Biopython is faster, often substantially. Use it.

Measured on an M-series Mac, hg38, one core, warm caches. Reproducible via the scripts in `benchmarks/`.

### Operations SeqMat is built for

| Operation                                          | Time     |
| -------------------------------------------------- | -------: |
| `gene_names_at_position(chrm, pos)`                | 2.5 µs   |
| KRAS mature mRNA assembly (splice + translate)     | 0.2 ms   |
| 1,000-SNP batch on 45 kb sequence (with history)   | 0.5 ms   |
| `Gene.from_file("KRAS")` (SQLite + unpickle)       | 24 ms    |
| `Gene.from_position(chrm, pos)` end-to-end         | 24 ms    |

### Position → gene, head-to-head

Same 63,241 hg38 gene intervals, same 10,000 random point queries. PyRanges is the natural comparison; both libraries are good at this, in different idioms. Reproduce with `python benchmarks/bench_position_lookup.py`.

Per-query (one coordinate, one answer, in a loop — `Gene.from_position`'s pattern):

| Implementation                          | Per query |
| --------------------------------------- | --------: |
| SeqMat (locator)                        | 2.5 µs    |
| Python `dict` + `bisect`                | 21 µs     |
| pandas (`groupby` chrm + boolean mask)  | 79 µs     |
| PyRanges constructed per call           | 2.0 ms    |

Batched (all queries handed in at once — PyRanges' native idiom):

| Implementation                | Per query |
| ----------------------------- | --------: |
| PyRanges `.join`              | 2.07 µs   |
| SeqMat (serial loop)          | 2.59 µs   |

Each library wins its native workload. Constructing a `PyRanges` per call is slow (and not what its authors intend) — listed only because it's a common mistake.

### Where Biopython is faster

For completeness, here are the operations where SeqMat is the wrong choice. The gap is the cost of the structured-array data model, which carries `(nt, ref, index, mut_type, valid)` through every transform.

| Workload                       | Biopython / native              | SeqMat  |
| ------------------------------ | ------------------------------: | ------: |
| 1,000 SNPs (no history needed) | `MutableSeq`: 0.21 ms           | 0.5 ms (with history) |
| 10-intron splice (50 kb)       | `str.join`: 0.002 ms            | 0.8 ms (coordinate-preserving) |
| Reverse-complement (1 Mb)      | `Seq.reverse_complement`: 0.57 ms | 9.4 ms |

If your hot path is one of these and you don't need coordinates or history attached to the result, use `Bio.Seq` or `bytearray`. The two libraries compose fine.

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

## Roadmap

[ROADMAP.md](ROADMAP.md) — what's coming, what isn't, and why. The headline themes: VCF + per-patient `Patient` / `Cohort` objects with lazy region-scoped variant application, reference-build-aware coordinate conversion, and splicing-impact integration.

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
