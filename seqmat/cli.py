#!/usr/bin/env python3
"""SeqMat CLI: setup, organisms, summary, biotypes, count, list, search, info, build-lmdb, test."""

import argparse
import sys

from .config import (
    get_available_organisms,
    get_default_organism,
    get_data_dir,
    get_organism_info as get_organism_config_info,
)
from .lmdb_store import build_lmdb
from .utils import (
    setup_genomics_data,
    print_data_summary,
    list_available_organisms,
    list_supported_organisms,
    list_gene_biotypes,
    count_genes,
    get_gene_list,
    search_genes,
    get_organism_info,
    test_installation,
)


def cmd_setup(args):
    """Setup genomics data for an organism"""
    try:
        setup_genomics_data(
            basepath=args.path,
            organism=args.organism,
            force=args.force,
            pickup=args.pickup
        )
        print(f"✅ Successfully set up {args.organism} data in {args.path}")
    except Exception as e:
        print(f"❌ Error setting up data: {e}")
        sys.exit(1)


def cmd_list_organisms(args):
    """List supported/configured organisms."""
    print("🌍 Organism Support Status:")
    print("-" * 30)
    supported = list_supported_organisms()
    configured = list_available_organisms()
    organism_names = {}
    for org in set(supported + configured):
        try:
            info = get_organism_config_info(org)
            organism_names[org] = info.get("name", org)
        except Exception:
            organism_names[org] = org
    
    for org in supported:
        name = organism_names.get(org, org)
        status = "✅ Configured" if org in configured else "❌ Not configured"
        print(f"{org}: {name} - {status}")
    
    if not configured:
        print("\nTo set up data, run:")
        print("  seqmat-setup --path /your/data/path --organism hg38")


def cmd_summary(args):
    """Print data summary"""
    print_data_summary()


def cmd_biotypes(args):
    """List gene biotypes for an organism"""
    if not args.organism:
        print("❌ Please specify an organism with --organism")
        sys.exit(1)
    
    biotypes = list_gene_biotypes(args.organism)
    
    if not biotypes:
        print(f"❌ No data found for organism '{args.organism}'")
        print("Available organisms:", ", ".join(list_available_organisms()))
        sys.exit(1)
    
    print(f"📊 Gene biotypes in {args.organism}:")
    print("-" * 30)
    counts = count_genes(args.organism)
    
    for biotype in biotypes:
        count = counts.get(biotype, 0)
        print(f"{biotype}: {count:,} genes")


def cmd_count(args):
    """Count genes for an organism/biotype"""
    if not args.organism:
        print("❌ Please specify an organism with --organism")
        sys.exit(1)
    
    counts = count_genes(args.organism, args.biotype)
    
    if not counts:
        print(f"❌ No data found for organism '{args.organism}'")
        sys.exit(1)
    
    if args.biotype:
        count = counts.get(args.biotype, 0)
        print(f"📊 {args.organism} {args.biotype}: {count:,} genes")
    else:
        print(f"📊 Gene counts for {args.organism}:")
        print("-" * 30)
        total = 0
        for biotype, count in sorted(counts.items()):
            print(f"{biotype}: {count:,} genes")
            total += count
        print("-" * 30)
        print(f"Total: {total:,} genes")


def cmd_list_genes(args):
    """List genes for an organism/biotype"""
    if not args.organism or not args.biotype:
        print("❌ Please specify both --organism and --biotype")
        sys.exit(1)
    
    genes = get_gene_list(args.organism, args.biotype, limit=args.limit)
    
    if not genes:
        print(f"❌ No genes found for {args.organism} {args.biotype}")
        sys.exit(1)
    
    print(f"📋 {args.organism} {args.biotype} genes ({len(genes)} shown):")
    print("-" * 50)
    
    for i, gene in enumerate(genes, 1):
        print(f"{i:4d}. {gene}")
    
    if args.limit and len(genes) == args.limit:
        total_count = count_genes(args.organism, args.biotype)
        total = total_count.get(args.biotype, 0)
        print(f"\n(Showing first {args.limit} of {total:,} total genes)")


def cmd_search(args):
    """Search for genes by name pattern"""
    if not args.organism or not args.query:
        print("❌ Please specify both --organism and --query")
        sys.exit(1)
    
    results = search_genes(
        organism=args.organism,
        query=args.query,
        biotype=args.biotype,
        limit=args.limit
    )
    
    if not results:
        print(f"❌ No genes found matching '{args.query}' in {args.organism}")
        sys.exit(1)
    
    print(f"🔍 Search results for '{args.query}' in {args.organism}:")
    print("-" * 50)
    
    for i, result in enumerate(results, 1):
        gene_id = result.get('gene_id', '')
        if gene_id:
            print(f"{i:2d}. {result['gene_name']} ({gene_id}) - {result['biotype']}")
        else:
            print(f"{i:2d}. {result['gene_name']} ({result['biotype']})")
    
    if len(results) == args.limit:
        print(f"\n(Showing first {args.limit} results)")


def cmd_build_lmdb(args):
    """Build LMDB database from gene pickle files (requires pip install seqmat[lmdb])"""
    try:
        output = build_lmdb(
            annotations_dir=args.annotations_dir,
            output_path=args.output,
            organism=args.organism,
        )
        if args.set_config:
            from .config import load_config, save_config
            config = load_config()
            org = args.organism or get_default_organism()
            if org in config and isinstance(config[org], dict):
                config[org]["gene_lmdb_path"] = output
            else:
                config["gene_lmdb_path"] = output
            save_config(config)
            print(f"Config updated: gene_lmdb_path = {output}")
    except ImportError as e:
        print(f"LMDB not installed: {e}")
        print("Install with: pip install seqmat[lmdb]")
        sys.exit(1)
    except Exception as e:
        print(f"Error building LMDB: {e}")
        sys.exit(1)


def cmd_info(args):
    """Show detailed information about an organism"""
    if not args.organism:
        print("❌ Please specify an organism with --organism")
        sys.exit(1)
    
    info = get_organism_info(args.organism)
    
    if "error" in info:
        print(f"❌ {info['error']}")
        sys.exit(1)
    
    print(f"ℹ️  Detailed information for {args.organism}:")
    print("=" * 40)
    data_avail = info.get("data_available", {})
    
    if "gene_counts" in data_avail:
        print("📊 Gene Data:")
        total_genes = 0
        for biotype, count in sorted(data_avail["gene_counts"].items()):
            print(f"  {biotype}: {count:,} genes")
            total_genes += count
        print(f"  Total: {total_genes:,} genes")
        print()
    
    if "chromosomes" in data_avail:
        chroms = data_avail["chromosomes"]
        print(f"🧬 Chromosome Data: {len(chroms)} chromosomes")
        print(f"  Available: {', '.join(sorted(chroms))}")
        print()
    
    print("📁 Data Paths:")
    for path_name, path_value in info["paths"].items():
        from pathlib import Path
        exists = "✅" if Path(path_value).exists() else "❌"
        print(f"  {path_name}: {exists} {path_value}")


def cmd_test(args):
    """Run SeqMat installation tests."""
    results = test_installation(args.organism, verbose=not args.quiet)
    sys.exit(1 if results["tests_failed"] > 0 else 0)


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(prog="seqmat", description="SeqMat genomics data management CLI")
    subparsers = parser.add_subparsers(dest="command", help="Commands")

    setup_parser = subparsers.add_parser("setup", help="Set up genomics data")
    default_data_dir = str(get_data_dir())
    setup_parser.add_argument("--path", default=default_data_dir,
                              help=f"Base path for data storage (default: {default_data_dir})")
    available_organisms = get_available_organisms()
    default_organism = get_default_organism()
    setup_parser.add_argument("--organism", default=default_organism, choices=available_organisms,
                             help=f"Organism to set up (default: {default_organism})")
    setup_parser.add_argument("--force", action="store_true", help="Force overwrite existing data")
    setup_parser.add_argument("--pickup", action="store_true", help="Resume interrupted setup, reuse existing downloaded files")
    setup_parser.set_defaults(func=cmd_setup)

    organisms_parser = subparsers.add_parser("organisms", help="List supported/configured organisms")
    organisms_parser.set_defaults(func=cmd_list_organisms)

    summary_parser = subparsers.add_parser("summary", help="Show data summary")
    summary_parser.set_defaults(func=cmd_summary)

    biotypes_parser = subparsers.add_parser("biotypes", help="List gene biotypes")
    biotypes_parser.add_argument("--organism", help="Organism to query")
    biotypes_parser.set_defaults(func=cmd_biotypes)

    count_parser = subparsers.add_parser("count", help="Count genes")
    count_parser.add_argument("--organism", help="Organism to query")
    count_parser.add_argument("--biotype", help="Specific biotype to count")
    count_parser.set_defaults(func=cmd_count)

    list_parser = subparsers.add_parser("list", help="List genes")
    list_parser.add_argument("--organism", help="Organism to query")
    list_parser.add_argument("--biotype", help="Gene biotype")
    list_parser.add_argument("--limit", type=int, default=50, help="Maximum genes to show")
    list_parser.set_defaults(func=cmd_list_genes)

    search_parser = subparsers.add_parser("search", help="Search genes by name")
    search_parser.add_argument("--organism", help="Organism to search")
    search_parser.add_argument("--query", help="Gene name pattern to search")
    search_parser.add_argument("--biotype", help="Filter by biotype")
    search_parser.add_argument("--limit", type=int, default=20, help="Maximum results")
    search_parser.set_defaults(func=cmd_search)

    info_parser = subparsers.add_parser("info", help="Show organism information")
    info_parser.add_argument("--organism", help="Organism to query")
    info_parser.set_defaults(func=cmd_info)

    lmdb_parser = subparsers.add_parser("build-lmdb", help="Build LMDB from gene pickle files (pip install seqmat[lmdb])")
    lmdb_parser.add_argument("--organism", default=None, help="Organism (default from config)")
    lmdb_parser.add_argument("--annotations-dir", default=None, help="Annotations dir (default from organism config)")
    lmdb_parser.add_argument("--output", default=None, help="Output LMDB path (default: <annotations>/genes.lmdb)")
    lmdb_parser.add_argument("--set-config", action="store_true", help="Update config to use new LMDB")
    lmdb_parser.set_defaults(func=cmd_build_lmdb)

    test_parser = subparsers.add_parser("test", help="Test SeqMat installation and data setup")
    test_parser.add_argument("--organism", help="Organism to test (uses default if not specified)")
    test_parser.add_argument("--quiet", action="store_true", help="Suppress detailed output")
    test_parser.set_defaults(func=cmd_test)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()