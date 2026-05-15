"""Honest benchmark of core sequence operations vs Biopython and naive baselines.

Three workloads — the operations SeqMat is built around:

1. **SNP application** with mutation tracking — apply 1,000 SNPs to a 50 kb sequence.
2. **Splicing** — excise 10 intron-sized regions from a 50 kb sequence.
3. **Reverse complement** — flip a 1 Mb sequence.

Each comparison is on the same underlying data. ``construct`` cost is reported
separately so that the "apply" line is what the call site would actually pay.

Reproduce:
    python benchmarks/bench_sequence_ops.py
    python benchmarks/bench_sequence_ops.py --runs 10
"""
from __future__ import annotations

import argparse
import random
import time
from typing import Callable

from Bio.Seq import MutableSeq, Seq
from seqmat import SeqMat


def best_of(fn: Callable[[], object], runs: int) -> float:
    best = float("inf")
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def fmt(seconds: float) -> str:
    if seconds < 1e-3:
        return f"{seconds * 1e6:7.1f} us"
    if seconds < 1.0:
        return f"{seconds * 1e3:7.2f} ms"
    return f"{seconds:7.2f}  s"


def header(title: str) -> None:
    print(f"\n{title}")
    print("-" * len(title))


def row(label: str, t: float, baseline: float | None) -> None:
    rel = f"{t / baseline:6.1f}x" if baseline else "1×"
    print(f"  {label:<46} {fmt(t):>11}   {rel:>8}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--runs", type=int, default=5)
    args = p.parse_args()
    runs = args.runs

    rng = random.Random(0)
    seq_50k = "ATCG" * 12_500
    seq_1m = "".join(rng.choices("ACGT", k=1_000_000))

    # ----- 1. SNP application -----
    N = 1_000
    zero_positions = sorted(random.Random(0).sample(range(len(seq_50k)), N))
    new_bases = [rng.choice("ACGT") for _ in range(N)]
    one_positions = [p + 1 for p in zero_positions]            # SeqMat is 1-indexed by default

    template = SeqMat(seq_50k)
    seq_str = template.seq
    muts = list(zip(one_positions, [seq_str[p] for p in zero_positions], new_bases))

    def sm_apply():
        s = SeqMat(seq_50k)
        s.apply_mutations(muts, permissive_ref=True)

    def bio_apply():
        s = MutableSeq(seq_50k)
        for p, b in zip(zero_positions, new_bases):
            s[p] = b

    def bytearray_apply():
        s = bytearray(seq_50k.encode())
        for p, b in zip(zero_positions, new_bases):
            s[p] = ord(b)

    t_construct_50k = best_of(lambda: SeqMat(seq_50k), runs)

    header(f"1) Apply {N:,} SNPs to a 50 kb sequence")
    t_sm  = best_of(sm_apply, runs) - t_construct_50k
    t_bio = best_of(bio_apply, runs)
    t_ba  = best_of(bytearray_apply, runs)
    row("Biopython MutableSeq (no tracking)", t_bio, t_bio)
    row("bytearray (no tracking)", t_ba, t_bio)
    row("SeqMat (apply only, with mutation history)", t_sm, t_bio)
    print(f"    (SeqMat construct cost subtracted: {fmt(t_construct_50k)})")

    # ----- 2. Splicing -----
    introns = [(2_000, 5_000), (8_000, 11_000), (15_000, 17_500),
               (20_000, 22_500), (25_500, 28_000), (31_000, 33_500),
               (36_000, 38_500), (41_000, 43_500), (45_500, 47_000),
               (48_000, 49_000)]

    def sm_splice():
        s = SeqMat(seq_50k)
        s.remove_regions(introns)

    def py_splice():
        keep = []
        last = 0
        for a, b in introns:
            keep.append(seq_50k[last:a])
            last = b
        keep.append(seq_50k[last:])
        return "".join(keep)

    header("2) Splice 10 introns out of a 50 kb sequence")
    t_sm = best_of(sm_splice, runs) - t_construct_50k
    t_py = best_of(py_splice, runs)
    row("Python str slice + join", t_py, t_py)
    row("SeqMat.remove_regions (coordinate-tracked)", t_sm, t_py)

    # ----- 3. Reverse complement -----
    def sm_rc():
        s = SeqMat(seq_1m)
        s.reverse_complement()

    def bio_rc():
        return str(Seq(seq_1m).reverse_complement())

    _RC = bytes.maketrans(b"ACGT", b"TGCA")

    def bytes_rc():
        return seq_1m.encode().translate(_RC)[::-1].decode()

    t_construct_1m = best_of(lambda: SeqMat(seq_1m), runs)

    header("3) Reverse-complement a 1 Mb sequence")
    t_bio = best_of(bio_rc, runs)
    t_by  = best_of(bytes_rc, runs)
    t_sm  = best_of(sm_rc, runs) - t_construct_1m
    row("Biopython Seq.reverse_complement", t_bio, t_bio)
    row("bytes.translate + [::-1]", t_by, t_bio)
    row("SeqMat.reverse_complement (coordinate-tracked)", t_sm, t_bio)
    print(f"    (SeqMat construct cost subtracted: {fmt(t_construct_1m)})")

    print("\nNote: Biopython / bytearray / str baselines don't track mutation history,")
    print("don't preserve genomic coordinates through indels, and don't keep the")
    print("reference alongside the current sequence — features only SeqMat provides.")


if __name__ == "__main__":
    main()
