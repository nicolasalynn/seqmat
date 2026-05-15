"""Compare chromosome+position → gene lookup against pyranges and pandas.

Honest, apples-to-apples:
- All implementations index the same ``(chrm, start, end, name)`` rows from ``genes.db``.
- Setup time (one-time index build) is timed separately from steady-state queries.
- Steady-state runs 10,000 random point queries drawn from the gene set so every
  query hits a real interval.
- Each timing is the *minimum* of 5 runs to suppress GC/JIT noise.

Usage:
    python benchmarks/bench_position_lookup.py
    python benchmarks/bench_position_lookup.py --organism hg38 --queries 10000 --runs 5
"""
from __future__ import annotations

import argparse
import pickle
import random
import sqlite3
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Tuple

import numpy as np

from seqmat import gene_names_at_position
from seqmat.config import get_default_organism
from seqmat.locator import _get_index
from seqmat.sqlite_store import get_genes_db_path


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

@dataclass
class GeneRow:
    chrm: str
    start: int
    end: int
    name: str


def load_gene_rows(organism: str) -> List[GeneRow]:
    """Extract (chrm, start, end, name) from genes.db."""
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
    """Sample N random (chrm, pos) queries from real gene midpoints."""
    rng = random.Random(seed)
    sampled = [rng.choice(rows) for _ in range(n)]
    return [(r.chrm, (r.start + r.end) // 2) for r in sampled]


# ---------------------------------------------------------------------------
# Implementations
# ---------------------------------------------------------------------------

def build_seqmat(organism: str) -> Callable[[str, int], list]:
    """Warm SeqMat's in-memory + sidecar index, return a query closure."""
    _get_index(organism)  # warm
    def query(chrm: str, pos: int) -> list:
        return gene_names_at_position(chrm, pos, organism=organism)
    return query


def build_pyranges(rows: List[GeneRow]):
    """Build a PyRanges from the same row set; return a query closure."""
    import pandas as pd
    import pyranges as pr

    df = pd.DataFrame({
        "Chromosome": [r.chrm for r in rows],
        "Start":      [r.start for r in rows],
        "End":        [r.end for r in rows],
        "Name":       [r.name for r in rows],
    })
    gr = pr.PyRanges(df)

    def query(chrm: str, pos: int) -> list:
        # PyRanges expects "chr"-prefixed names if the index was built that way.
        # Our rows are stored without prefix; mirror that on the query side.
        q = pr.PyRanges(pd.DataFrame({"Chromosome": [chrm], "Start": [pos], "End": [pos + 1]}))
        hit = gr.join(q)
        return hit.Name.tolist() if len(hit) else []

    return query


def build_pandas(rows: List[GeneRow]):
    """Naive pandas baseline: filter on chrm + boolean mask on the full DataFrame."""
    import pandas as pd

    df = pd.DataFrame({
        "chrm":  [r.chrm for r in rows],
        "start": [r.start for r in rows],
        "end":   [r.end for r in rows],
        "name":  [r.name for r in rows],
    })

    def query(chrm: str, pos: int) -> list:
        hit = df[(df.chrm == chrm) & (df.start <= pos) & (df.end >= pos)]
        return hit.name.tolist()

    return query


def build_dict_bisect(rows: List[GeneRow]):
    """Per-chromosome sorted starts + bisect baseline — what someone might write by hand."""
    import bisect
    from collections import defaultdict

    per_chrm = defaultdict(list)
    for r in rows:
        per_chrm[r.chrm].append((r.start, r.end, r.name))
    for c in per_chrm:
        per_chrm[c].sort(key=lambda t: t[0])
    starts_by_chrm = {c: [t[0] for t in v] for c, v in per_chrm.items()}

    def query(chrm: str, pos: int) -> list:
        entries = per_chrm.get(chrm)
        if not entries:
            return []
        starts = starts_by_chrm[chrm]
        hi = bisect.bisect_right(starts, pos)
        return [name for (s, e, name) in entries[:hi] if e >= pos]

    return query


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def best_of(fn: Callable[[], None], runs: int) -> float:
    """Return the minimum wall-clock seconds over ``runs`` invocations."""
    best = float("inf")
    for _ in range(runs):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


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

    implementations: List[Tuple[str, Callable]] = []

    print("\nBuilding indexes (one-time setup cost):")

    t0 = time.perf_counter()
    seqmat_q = build_seqmat(organism)
    print(f"  SeqMat (cached sidecar):  {(time.perf_counter() - t0) * 1e3:.1f} ms")
    implementations.append(("SeqMat", seqmat_q))

    t0 = time.perf_counter()
    pyranges_q = build_pyranges(rows)
    print(f"  PyRanges (PyRanges):      {(time.perf_counter() - t0) * 1e3:.1f} ms")
    implementations.append(("PyRanges 0.1.x", pyranges_q))

    t0 = time.perf_counter()
    pandas_q = build_pandas(rows)
    print(f"  pandas DataFrame:         {(time.perf_counter() - t0) * 1e3:.1f} ms")
    implementations.append(("pandas (boolean mask)", pandas_q))

    t0 = time.perf_counter()
    bisect_q = build_dict_bisect(rows)
    print(f"  dict + bisect (handroll): {(time.perf_counter() - t0) * 1e3:.1f} ms")
    implementations.append(("Python dict + bisect", bisect_q))

    print("\nSteady-state query benchmark (lower = better):")
    print(f"{'Implementation':<28} {'Total (s)':>10}  {'Per query':>11}  {'Relative':>10}")
    print("-" * 64)

    # Warm each impl once (avoid cold-cache surprises)
    for _, q in implementations:
        for chrm, pos in queries[:50]:
            q(chrm, pos)

    seqmat_per = None
    results = []
    for label, q in implementations:
        # PyRanges per-query is dramatically slower; cap its query count so the run finishes in seconds.
        # We still report PER-QUERY latency, which is the apples-to-apples number.
        if label.startswith("PyRanges"):
            sample = queries[:max(50, args.queries // 100)]
        else:
            sample = queries

        def run(qs=sample, fn=q):
            for chrm, pos in qs:
                fn(chrm, pos)

        total = best_of(run, runs=args.runs)
        per_us = total / len(sample) * 1e6
        results.append((label, total, per_us, len(sample)))
        if label == "SeqMat":
            seqmat_per = per_us

    for label, total, per_us, n in results:
        rel = f"{per_us / seqmat_per:>7.1f}x" if seqmat_per else "—"
        if per_us > 1_000:
            per_str = f"{per_us / 1000:6.2f} ms"
        else:
            per_str = f"{per_us:6.2f} us"
        print(f"{label:<28} {total:>10.3f}  {per_str:>11}  {rel:>10}    (n={n:,})")


if __name__ == "__main__":
    main()
