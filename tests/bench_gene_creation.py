"""Benchmark: Gene + Transcript creation throughput.

Measures the time to load a gene, pull pre-mRNA, generate mature mRNA,
and translate to protein for clinically important genes.

Usage:
    python tests/bench_gene_creation.py
    python tests/bench_gene_creation.py --genes TP53 BRCA1 --iterations 50
"""
import argparse
import statistics
import time
from typing import Dict, List, Tuple

from seqmat import Gene

DEFAULT_GENES = ["TP53", "BRCA1", "BRCA2", "EGFR", "KRAS", "BRAF", "PIK3CA", "APC", "PTEN", "MYC"]
DEFAULT_ITERATIONS = 100


def bench_full_pipeline(gene_name: str, iterations: int) -> Dict[str, List[float]]:
    """Time each stage of the gene→transcript→protein pipeline."""
    timings: Dict[str, List[float]] = {
        "from_file": [],
        "generate_pre_mrna": [],
        "generate_mature_mrna": [],
        "generate_protein": [],
        "total": [],
    }

    # Warm-up (1 call to prime caches / FASTA index)
    g = Gene.from_file(gene_name)
    if g is None:
        print(f"  SKIP {gene_name}: not found in data")
        return timings
    t = g.transcript()
    if t is None:
        print(f"  SKIP {gene_name}: no primary transcript")
        return timings
    t.generate_mature_mrna()
    t.generate_protein()

    for _ in range(iterations):
        t0 = time.perf_counter()

        t1 = time.perf_counter()
        g = Gene.from_file(gene_name)
        t2 = time.perf_counter()
        timings["from_file"].append(t2 - t1)

        # Transcript.__init__ already calls generate_pre_mrna
        t3 = time.perf_counter()
        tx = g.transcript()
        t4 = time.perf_counter()
        timings["generate_pre_mrna"].append(t4 - t3)

        t5 = time.perf_counter()
        tx.generate_mature_mrna()
        t6 = time.perf_counter()
        timings["generate_mature_mrna"].append(t6 - t5)

        t7 = time.perf_counter()
        tx.generate_protein()
        t8 = time.perf_counter()
        timings["generate_protein"].append(t8 - t7)

        timings["total"].append(t8 - t0)

    return timings


def fmt_ms(values: List[float]) -> Tuple[str, str, str]:
    """Return (mean, median, stdev) formatted in ms."""
    if not values:
        return ("N/A", "N/A", "N/A")
    ms = [v * 1000 for v in values]
    mean = statistics.mean(ms)
    med = statistics.median(ms)
    sd = statistics.stdev(ms) if len(ms) > 1 else 0.0
    return (f"{mean:.2f}", f"{med:.2f}", f"{sd:.2f}")


def main():
    parser = argparse.ArgumentParser(description="Benchmark seqmat gene creation pipeline")
    parser.add_argument("--genes", nargs="+", default=DEFAULT_GENES, help="Gene names to benchmark")
    parser.add_argument("--iterations", "-n", type=int, default=DEFAULT_ITERATIONS, help="Iterations per gene")
    args = parser.parse_args()

    col_w = 14
    stages = ["from_file", "generate_pre_mrna", "generate_mature_mrna", "generate_protein", "total"]
    header_stages = ["load", "pre_mrna", "mature_mrna", "protein", "TOTAL"]

    # Header
    print(f"\nBenchmark: {args.iterations} iterations per gene\n")
    print(f"{'Gene':<10}", end="")
    for h in header_stages:
        print(f"  {h + ' (ms)':>{col_w}}", end="")
    print()
    print("-" * (10 + (col_w + 2) * len(stages)))

    all_totals: List[float] = []

    for gene_name in args.genes:
        timings = bench_full_pipeline(gene_name, args.iterations)
        if not timings["total"]:
            continue
        all_totals.extend(timings["total"])
        print(f"{gene_name:<10}", end="")
        for stage in stages:
            mean, _, _ = fmt_ms(timings[stage])
            print(f"  {mean:>{col_w}}", end="")
        print()

    # Summary
    if all_totals:
        print("-" * (10 + (col_w + 2) * len(stages)))
        ms_totals = [v * 1000 for v in all_totals]
        n_genes = len([g for g in args.genes if any(True for _ in [])] or args.genes)
        total_genes = len(args.genes)
        print(f"\nOverall across {total_genes} genes x {args.iterations} iterations:")
        print(f"  Mean total:   {statistics.mean(ms_totals):.2f} ms")
        print(f"  Median total: {statistics.median(ms_totals):.2f} ms")
        print(f"  Stdev:        {statistics.stdev(ms_totals):.2f} ms")
        print(f"  Min:          {min(ms_totals):.2f} ms")
        print(f"  Max:          {max(ms_totals):.2f} ms")
        print(f"  Throughput:   {1000 / statistics.mean(ms_totals):.1f} genes/sec")


if __name__ == "__main__":
    main()
