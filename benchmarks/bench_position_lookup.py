"""Honest comparison of chromosome+position → gene lookup.

Two workloads, two tables:

1. **Per-query (serial loop)** — call ``query(chrm, pos)`` 10,000 times in a row,
   exactly the way ``Gene.from_position`` is meant to be used. Reflects "I just
   parsed one VCF line, what gene is here?"
2. **Batched** — hand the implementation all 10,000 (chrm, pos) pairs at once,
   measure the total. Reflects "I have a full VCF, annotate everything."

Why both: PyRanges is built for batch interval-set algebra; calling it
per-query is the wrong tool. SeqMat's locator is built for per-query
lookups and doesn't have a vectorized batch API today. Each library wins
its native workload — the table tells the truth.

All implementations index the same ``(chrm, start, end, name)`` rows
extracted from ``genes.db``. Setup time (one-time index build) is timed
separately. Each timing is the minimum of N runs to suppress GC noise.

Usage:
    python benchmarks/bench_position_lookup.py
    python benchmarks/bench_position_lookup.py --queries 10000 --runs 5
"""
from __future__ import annotations

import argparse
import bisect
import pickle
import random
import sqlite3
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Tuple

import numpy as np

from seqmat import gene_names_at_position
from seqmat.config import get_default_organism
from seqmat.locator import _get_index
from seqmat.sqlite_store import get_genes_db_path


@dataclass
class GeneRow:
    chrm: str
    start: int
    end: int
    name: str


def load_gene_rows(organism: str) -> List[GeneRow]:
    db_path = get_genes_db_path(organism)
    if not db_path:
        raise RuntimeError(f"genes.db not found for organism {organism!r}")
    rows: List[GeneRow] = []
    with sqlite3.connect(db_path) as conn:
        for gene_name, gene_id, blob in conn.execute("SELECT gene_name, gene_id, data FROM genes"):
            data = pickle.loads(blob)
            chrm = data.get("chrm")
            start = data.get("gene_start")
            end = data.get("gene_end")
            if chrm is None or start is None or end is None:
                continue
            if start > end:
                start, end = end, start
            rows.append(GeneRow(str(chrm), int(start), int(end), gene_name or gene_id))
    return rows


def sample_queries(rows: List[GeneRow], n: int, seed: int = 0) -> List[Tuple[str, int]]:
    rng = random.Random(seed)
    sampled = [rng.choice(rows) for _ in range(n)]
    return [(r.chrm, (r.start + r.end) // 2) for r in sampled]


def best_of(fn: Callable[[], object], runs: int) -> float:
    best = float("inf")
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


# ---------------------------------------------------------------------------
# Indexes
# ---------------------------------------------------------------------------

def build_seqmat_index(organism: str):
    _get_index(organism)
    return organism


def build_pyranges_index(rows: List[GeneRow]):
    import pandas as pd
    import pyranges as pr
    df = pd.DataFrame({
        "Chromosome": [r.chrm for r in rows],
        "Start":      [r.start for r in rows],
        "End":        [r.end for r in rows],
        "Name":       [r.name for r in rows],
    })
    return pr.PyRanges(df)


def build_pandas_index(rows: List[GeneRow]):
    import pandas as pd
    df = pd.DataFrame({
        "chrm":  [r.chrm for r in rows],
        "start": np.array([r.start for r in rows], dtype=np.int64),
        "end":   np.array([r.end for r in rows], dtype=np.int64),
        "name":  [r.name for r in rows],
    })
    return {c: g for c, g in df.groupby("chrm")}


def build_dict_bisect_index(rows: List[GeneRow]):
    per_chrm = defaultdict(list)
    for r in rows:
        per_chrm[r.chrm].append((r.start, r.end, r.name))
    for c in per_chrm:
        per_chrm[c].sort(key=lambda t: t[0])
    starts_by_chrm = {c: [t[0] for t in v] for c, v in per_chrm.items()}
    return per_chrm, starts_by_chrm


# ---------------------------------------------------------------------------
# Per-query (serial loop) implementations
# ---------------------------------------------------------------------------

def serial_seqmat(organism, queries):
    for chrm, pos in queries:
        gene_names_at_position(chrm, pos, organism=organism)


def serial_dict_bisect(idx, queries):
    per_chrm, starts_by_chrm = idx
    for chrm, pos in queries:
        entries = per_chrm.get(chrm)
        if not entries:
            continue
        starts = starts_by_chrm[chrm]
        hi = bisect.bisect_right(starts, pos)
        [name for (s, e, name) in entries[:hi] if e >= pos]


def serial_pandas(by_chrm, queries):
    for chrm, pos in queries:
        g = by_chrm.get(chrm)
        if g is None:
            continue
        g.loc[(g.start <= pos) & (g.end >= pos), "name"].tolist()


def serial_pyranges(gr, queries):
    """Anti-pattern: constructing a PyRanges per call. Included only as a warning."""
    import pandas as pd
    import pyranges as pr
    for chrm, pos in queries:
        q = pr.PyRanges(pd.DataFrame({"Chromosome": [chrm], "Start": [pos], "End": [pos + 1]}))
        gr.join(q)


# ---------------------------------------------------------------------------
# Batched implementations
# ---------------------------------------------------------------------------

def batched_pyranges(gr, queries):
    """Construct one PyRanges of all queries, do one .join(). The natural PyRanges idiom."""
    import pandas as pd
    import pyranges as pr
    qdf = pd.DataFrame({
        "Chromosome": [q[0] for q in queries],
        "Start":      [q[1] for q in queries],
        "End":        [q[1] + 1 for q in queries],
    })
    qpr = pr.PyRanges(qdf)
    return gr.join(qpr)


def batched_seqmat(organism, queries):
    """SeqMat has no vectorized batch API today — this is the same serial loop."""
    for chrm, pos in queries:
        gene_names_at_position(chrm, pos, organism=organism)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--organism", default=None)
    p.add_argument("--queries", type=int, default=10_000)
    p.add_argument("--runs", type=int, default=5)
    args = p.parse_args()

    organism = args.organism or get_default_organism()
    print(f"Loading gene rows for {organism}...")
    rows = load_gene_rows(organism)
    print(f"  {len(rows):,} intervals across {len({r.chrm for r in rows})} chromosomes")

    print(f"Sampling {args.queries:,} random point queries (seeded)...")
    queries = sample_queries(rows, args.queries)

    print("\nBuilding indexes:")
    t0 = time.perf_counter()
    seqmat_idx = build_seqmat_index(organism)
    print(f"  SeqMat (cached sidecar):    {(time.perf_counter() - t0) * 1e3:6.1f} ms")
    t0 = time.perf_counter()
    pyranges_idx = build_pyranges_index(rows)
    print(f"  PyRanges:                   {(time.perf_counter() - t0) * 1e3:6.1f} ms")
    t0 = time.perf_counter()
    pandas_idx = build_pandas_index(rows)
    print(f"  pandas (groupby chrm):      {(time.perf_counter() - t0) * 1e3:6.1f} ms")
    t0 = time.perf_counter()
    dict_idx = build_dict_bisect_index(rows)
    print(f"  dict + bisect (handroll):   {(time.perf_counter() - t0) * 1e3:6.1f} ms")

    def fmt_per_query(per_us: float) -> str:
        return f"{per_us / 1000:6.2f} ms" if per_us > 1_000 else f"{per_us:6.2f} us"

    def print_row(label, total, per_us, rel_baseline, n):
        rel = f"{per_us / rel_baseline:>7.2f}x" if rel_baseline else "—"
        print(f"  {label:<32} {total*1e3:>7.1f} ms total   {fmt_per_query(per_us):>11}/query   {rel:>10}    (n={n:,})")

    # --- Per-query table ----------------------------------------------------
    print("\nPer-query (serial loop) — calling once per coordinate, the Gene.from_position pattern:")
    # Warm all
    for chrm, pos in queries[:50]:
        gene_names_at_position(chrm, pos, organism=organism)

    serial_results = []

    total = best_of(lambda: serial_seqmat(seqmat_idx, queries), args.runs)
    seqmat_per = total / len(queries) * 1e6
    serial_results.append(("SeqMat (locator)", total, seqmat_per, len(queries)))

    total = best_of(lambda: serial_dict_bisect(dict_idx, queries), args.runs)
    serial_results.append(("Python dict + bisect", total, total / len(queries) * 1e6, len(queries)))

    total = best_of(lambda: serial_pandas(pandas_idx, queries), max(1, args.runs // 2))
    serial_results.append(("pandas (groupby chrm)", total, total / len(queries) * 1e6, len(queries)))

    # PyRanges anti-pattern: only run a small sample so the script finishes
    sample = queries[:max(50, args.queries // 100)]
    total = best_of(lambda: serial_pyranges(pyranges_idx, sample), 1)
    serial_results.append(("PyRanges (constructed per call)", total, total / len(sample) * 1e6, len(sample)))

    for label, total, per_us, n in serial_results:
        print_row(label, total, per_us, seqmat_per, n)

    # --- Batched table ------------------------------------------------------
    print("\nBatched — all queries handed to the library at once (PyRanges' native idiom):")
    batch_results = []

    total = best_of(lambda: batched_pyranges(pyranges_idx, queries), args.runs)
    pyranges_batch_per = total / len(queries) * 1e6
    batch_results.append(("PyRanges (.join)", total, pyranges_batch_per, len(queries)))

    total = best_of(lambda: batched_seqmat(seqmat_idx, queries), args.runs)
    batch_results.append(("SeqMat (serial loop)", total, total / len(queries) * 1e6, len(queries)))

    for label, total, per_us, n in batch_results:
        print_row(label, total, per_us, pyranges_batch_per, n)


if __name__ == "__main__":
    main()
