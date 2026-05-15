# Changelog

All notable changes to **SeqMat** are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- `Gene.get(name)` — same lookup as `from_file` but raises typed exceptions
  (`GeneNotFoundError`, `OrganismNotConfiguredError`) instead of returning `None`.
  `Gene.from_file` is unchanged.
- Comparative benchmark (`benchmarks/bench_position_lookup.py`) vs PyRanges,
  pandas, and a handrolled `dict + bisect`. Numbers wired into the README.
- KRAS G12D end-to-end notebook (`examples/kras_g12d_analysis.ipynb`).
- `pdoc`-rendered API reference auto-deployed to GitHub Pages on every push.
- Repo logo, CI badge, Codecov badge.

## [1.5.0] — 2026-05-15

### Added
- Typed exception hierarchy in `seqmat.errors`: `SeqMatError` → `GeneNotFoundError`,
  `OrganismNotConfiguredError`, `DataUnavailableError`
  (`PrebuiltDataUnavailableError` kept as a backward-compatible alias).
- `CHANGELOG.md`, GitHub Actions CI (pytest on 3.10/3.11/3.12 + sdist/wheel build).
- `pyproject.toml` (PEP 621); `setup.py` removed.

### Changed
- `utils.py` (1,603 lines) split into focused modules: `_io.py`, `download.py`,
  `gtf.py`, `conservation.py`, `data_setup.py`, `discovery.py`. `utils.py` remains
  a re-export shim so every previous import path keeps working.
- Library code now uses `logging.getLogger(__name__)` instead of bare `print()`.
  CLI keeps user-facing prints. Callers control verbosity.

### Fixed
- `list_gene_biotypes`, `count_genes`, `get_gene_list`, `search_genes`,
  `available_genes`, and `seqmat summary` previously scanned a legacy `*.pkl`
  directory tree that doesn't exist on modern SQLite-backed installs and
  silently returned empty results. They now query `genes.db` directly and
  correctly report (e.g.) 42 biotypes / 120k genes across hg38 + mm39.

## [1.4.0] — 2026-05-15

### Added
- `Gene.from_position(chrm, pos, organism=None) → list[Gene]` for point and range
  position queries. Backed by a per-chromosome sorted NumPy index persisted as
  `gene_locations.npz` next to `genes.db`. ~17 µs per hot lookup.
- `gene_names_at_position()` and `build_location_index()` module-level helpers.
- Build pipeline now emits `gene_locations.npz` alongside `genes.db` automatically.

### Changed
- `sqlite_store.load_gene_from_sqlite` now matches on `gene_name` **or** `gene_id`,
  so Ensembl entries with empty gene symbols (lncRNAs/pseudogenes) resolve correctly
  from position queries.

## [1.3.1] — 2026-03-09

### Fixed
- File-handle leak in FASTA reader.
- Removed broken `Gene.__deepcopy__`.
- Dead code cleanup.

## [1.3.0] — 2026-03-09

### Changed
- Pre-mRNA generation is now lazy.
- Mature mRNA is assembled directly from exon spans instead of subtracting introns.

## [1.2.2] — 2026-03-09

### Fixed
- `clone()` on reverse-strand genes when `TIS > TTS`.

## [1.2.1] — 2026-03-09

### Changed
- Hot-path optimizations: `np.frombuffer`, LUT-based complement, `np.searchsorted`
  on sorted exon coordinates, ORF clone reuse.

## [1.2.0] — 2026-03-09

### Added
- FASTA caching layer for repeated range reads.
- Range-scoped reads (only the requested bytes are pulled).
- Copy-on-write `SeqMat.clone()`.

## [1.1.1] — 2026-02-28

### Changed
- Rebuilt wheel without `biopython` dependency.

## [1.1.0] — 2026-02-28

### Added
- `SeqMat.alignment()` — gapped reference/mutant sequence pair output.

### Removed
- `biopython` dependency. Standard codon translation table is now inlined.

## [1.0.0] — 2026-02-21

### Added
- First Production/Stable release.

## [0.x] — pre-1.0

Selected highlights from the pre-1.0 development history:

- **0.1.55** — Prebuilt S3 download default, config-less setup, build-from-sources path.
- **0.1.47** — `cons_vector` fallback, optional flanking, clearer ref-mismatch message.
- **0.1.45** — Pre-LMDB stable line; baseline for the production releases.

[Unreleased]: https://github.com/nicolasalynn/seqmat/compare/v1.5.0...HEAD
[1.5.0]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.5.0
[1.4.0]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.4.0
[1.3.1]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.3.1
[1.3.0]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.3.0
[1.2.2]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.2.2
[1.2.1]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.2.1
[1.2.0]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.2.0
[1.1.1]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.1.1
[1.1.0]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.1.0
[1.0.0]: https://github.com/nicolasalynn/seqmat/releases/tag/v1.0.0
