# SeqMat Roadmap

A living document for where the library is heading. Themes are independent; numbered phases at the bottom suggest an order.

## Vision

SeqMat occupies a narrow, well-defined slot: **the labeled-sequence layer for variant-effect and splicing analyses in Python.** Biopython holds the bases; pyensembl/pyranges hold the annotations; SeqMat holds the *coordinate-tracked, mutation-aware sequence* that sits between them and is currently always written by hand. The medium-term goal is to be the default answer to:

- *"I have a variant and a transcript — what's the protein effect?"*
- *"I have a VCF and a gene — give me the per-patient mature mRNA / protein."*
- *"Does this variant interact with a splice site?"*

These are common questions with no single Python-native answer today. They're what SeqMat's data model is shaped for.

---

## Theme 1 — VCF integration and the `Patient` / `Cohort` model

The single most-leveraged piece of work. Without VCF I/O, SeqMat isn't in the conversation. With it, the data model becomes immediately useful for clinical genomics workflows.

### Design sketch

```python
# Single-sample
patient = Patient.from_vcf("P001.vcf.gz")          # lazy: parse + index, do not apply
patient.region("12", 25_200_000, 25_300_000)       # SeqMat with variants in range applied
patient.gene("KRAS")                               # Gene whose pre_mRNA already has the variants
patient.transcript("KRAS").generate_protein()      # full pipeline, patient-specific

# Cohort
cohort = Cohort.from_vcf("cohort.vcf.gz")          # one VCF, N samples
cohort.samples                                      # ['P001', 'P002', ...]
cohort.gene("KRAS", sample="P002")                 # Gene for one patient
for sample_id, gene in cohort.iter_gene("KRAS"):
    ...                                            # streaming per-patient

# Vectorized cohort query
cohort.protein_effects("KRAS")                     # DataFrame: sample_id × { variants, protein_diff, ... }
```

### Why lazy

The data model is a `tabix`-indexed VCF on disk plus an in-memory index. When you ask for `patient.region("12", 25_200_000, 25_300_000)`, SeqMat:

1. Range-queries the VCF for variants in that span (single `tabix` call).
2. Builds the SeqMat for the region from FASTA.
3. Applies just those variants.

No global "apply all variants" pre-pass. A whole-exome VCF has millions of variants; pre-applying all of them is wasteful when most analyses touch a handful of genes.

### Implementation notes

- Depend on `cyvcf2` (fast pysam-backed VCF reader; multi-sample friendly).
- `Patient` is `(vcf_handle, sample_index, organism)`. Cheap to construct.
- `Cohort` wraps one multi-sample VCF — the common BCF/VCF format from joint genotyping. Streaming iteration matters here; a 1000-sample cohort × 30k genes shouldn't try to materialize.
- Variant-application caching per `(patient, region)` is optional; start without it and measure.

### Scope guard

Phased ploidy / phasing is a rabbit hole. v1 of `Patient` should default to "apply all variants where the sample carries the alt allele (hom or het), produce one sequence." Phasing-aware "give me each haplotype separately" can come later.

---

## Theme 2 — Killer differentiators

Two candidates that nobody else does cleanly. Pick one to lead on; both work.

### 2a. Reference-build-aware coordinate conversion

```python
variant = (12, 25_245_350, "C", "T")
SeqMat.lift("hg19", "hg38", variant)               # uses chain files, returns lifted variant
gene_hg19 = Gene.from_position("12", 25_245_350, organism="hg19")
gene_hg38 = gene_hg19.lift("hg38")                  # re-resolve from lifted coords
```

Reference-build mismatch is one of the most common sources of clinical-genomics bugs. Today: CrossMap + manual reasoning. SeqMat's coordinate model is exactly the right abstraction to make this safe by construction.

Dependencies: liftover chain files (UCSC), `pyliftover` or a direct chain-file reader.

### 2b. Splicing-impact integration

```python
tx = kras.transcript()
tx.predict_splice_effect((25_245_350, "C", "T"), tool="SpliceAI")
# returns: { score: 0.93, donor_loss: True, acceptor_gain: False, ... }
```

Plug SpliceAI / MaxEntScan / MMSplice in as a unified method on `Transcript`. The infrastructure for variant + transcript + position-relative-to-splice-site is already in SeqMat. Wrapping the predictors is the work.

This is the splicing-shop killer feature — it matches the repo's existing oncosplice references.

---

## Theme 3 — Interoperability

Removes the "do I have to migrate?" objection. Makes SeqMat composable with what people already use.

```python
gene = Gene.from_pyensembl(ensembl_release, "KRAS")
record: Bio.SeqRecord = transcript.to_seqrecord()
transcript = Transcript.from_seqrecord(record, ...)
seq_mat = SeqMat.from_seq(bio_seq, indices=...)
bio_seq = seq_mat.to_seq()
```

Cheap to implement, big legibility win. Lets SeqMat be *the layer above Biopython* rather than a fork.

---

## Theme 4 — The byte-storage architectural decision

The structured-array model is the source of both SeqMat's value (parallel columns) and its byte-op slowdown vs Biopython (27-byte records). Two paths:

### 4a. Delegate to Biopython internally

```
SeqMat
├─ self._mut:   MutableSeq          # Biopython, C-backed byte storage
├─ self._ref:   MutableSeq          # parallel reference
├─ self._index: np.ndarray[int64]   # genomic coordinates
├─ self._valid: np.ndarray[bool]    # splice mask
└─ self.mutations: list[dict]
```

`reverse_complement`, `complement`, `translate` delegate to `MutableSeq` — Biopython-speed for free.

Cost: hard dependency on Biopython; locality lost across columns (a per-row scan touches three objects instead of one). For most operations this is fine; for `_refresh_mutation_state`-style scans across `ref` + `nt` + `valid`, it's a regression.

### 4b. Cython/Rust the hot inner loops

Same data model, hot loops in compiled code. `reverse_complement_inplace`, `complement_inplace`, the SNP batch path, the structured-array reverse-copy. Likely a 5–10× speedup on byte ops while keeping the data model.

Cost: build wheel complexity, contributor friction.

**Recommendation:** Try 4a first. The pedagogical framing ("SeqMat is to Biopython what pandas is to numpy") is the strongest pitch for adoption, and the dependency is fair. Reach for 4b only if 4a doesn't close enough of the gap.

---

## Theme 5 — Documentation aimed at recipes

Right now the docs are API-oriented. The adoption-driving docs are recipe-oriented. Targets:

- "Annotate a VCF with per-variant protein effects"
- "Compare two reference builds for a variant"
- "Predict splice-site disruption for a coding variant"
- "Build a per-patient mature mRNA from a cohort VCF"
- "Find variants near splice sites across a gene panel"

Each recipe is ~50 lines and the README links to all of them. Recipe-driven docs win adoption because they answer the question the user actually has.

Side note: the existing `examples/` directory is the right home for these. They should be runnable, pre-executed notebooks like the KRAS one.

---

## Explicit non-goals

Things SeqMat should *not* try to do, even if asked:

- **General sequence I/O** beyond FASTA + VCF. No GenBank parser. No EMBL parser. Use `Bio.SeqIO`.
- **Pairwise or multiple-sequence alignment.** Use `Bio.Align`, `parasail`, `pyhmmer`.
- **Phylogenetics, structure (PDB), population genetics, restriction enzymes, motif analysis.** All Biopython modules with stable communities. Stay out.
- **BLAST / Entrez integration.** Use `Bio.Blast` / `Bio.Entrez`.
- **Beating `Bio.Seq.reverse_complement()` on long sequences.** The remaining gap doesn't matter for SeqMat-shaped workloads; pushing further is engineering for engineering's sake.

The narrative discipline: *if it's not a coordinate-aware or mutation-aware operation, it's not SeqMat's job.*

---

## Phasing

Suggested order of attack, with the work most likely to move the needle first.

**Phase A — VCF + Patient/Cohort (Theme 1)**
- `cyvcf2` dependency
- `Patient.from_vcf`, `patient.region`, `patient.gene`, `patient.transcript`
- `Cohort.from_vcf`, sample iteration
- A new notebook: *"From a VCF to a patient-specific protein"*

This alone makes SeqMat the only Python library that does this cleanly. It's also where the existing data model pays off most.

**Phase B — Interop (Theme 3)**
- `Gene.from_pyensembl`
- `to_seqrecord` / `from_seqrecord`
- `SeqMat.from_seq` / `SeqMat.to_seq`

Cheap, high-leverage, gives migration paths.

**Phase C — Killer feature: ref-build conversion (Theme 2a)**
- `SeqMat.lift` / `Gene.lift` over UCSC chain files
- Notebook: *"Does this hg19 variant mean the same thing in hg38?"*

**Phase D — Killer feature: splicing impact (Theme 2b)**
- `Transcript.predict_splice_effect`, wrapping SpliceAI initially
- Notebook in the splicing series

**Phase E — Architecture (Theme 4)**
- Decide between 4a and 4b based on real profiling against Phase-A workloads
- One release boundary later: SeqMat 2.0

**Phase F — Recipe docs (Theme 5)**
- One recipe per release until the catalog is solid

---

## Things I haven't decided

- Whether `Patient` should be a `Cohort` of size 1, or its own class. Probably its own class for ergonomics; `Cohort` becomes the multi-sample shape.
- Whether ref-build conversion belongs on `SeqMat` itself or on a free function. Probably both — `Gene.lift` is the headline API; `seqmat.lift(seq, "hg19", "hg38")` is the primitive.
- Whether to depend hard on `cyvcf2` (fast, but a C-extension) or accept a pure-Python fallback (`pysam.VariantFile`).

---

## Realistic outcomes (re-stated honestly)

- **Best case:** SeqMat is cited as the Python library for variant-effect-on-transcript work; ~1k active users; mentioned alongside Biopython in splicing/variant-effect papers.
- **Median case:** Useful niche tool for the author's group and a handful of others; quiet steady adoption.
- **Worst case:** A bigger lab releases something similar with funding; absorbs the niche.

The work that pushes from median toward best is **Theme 1 first, then 2a or 2b**. Two to three months of focused work, end to end.
