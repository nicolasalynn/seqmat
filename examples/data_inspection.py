#!/usr/bin/env python3
"""
Data Inspection and Management Examples

This script demonstrates the new data inspection utilities in SeqMat.
"""

from seqmat import (
    list_supported_organisms,
    list_available_organisms, 
    get_organism_info,
    list_gene_biotypes,
    count_genes,
    get_gene_list,
    search_genes,
    print_data_summary,
    data_summary
)

def basic_organism_info():
    """Demonstrate basic organism information retrieval"""
    print("=== Basic Organism Information ===")
    
    # Check what organisms are supported
    supported = list_supported_organisms()
    print(f"Supported organisms: {supported}")
    
    # Check what's actually configured
    available = list_available_organisms()
    print(f"Configured organisms: {available}")
    
    if not available:
        print("\n⚠️  No organisms configured yet!")
        print("Run: seqmat setup --path /your/data/path --organism hg38")
        return
    
    print()

def explore_organism_data():
    """Explore data for a specific organism"""
    print("=== Exploring Organism Data ===")
    
    available = list_available_organisms()
    if not available:
        print("No organisms configured. Run setup first.")
        return
    
    organism = available[0]  # Use first available organism
    print(f"Exploring data for: {organism}")
    
    # Get detailed organism information
    info = get_organism_info(organism)
    
    if "error" in info:
        print(f"Error: {info['error']}")
        return
    
    print(f"✅ Organism: {organism}")
    print(f"✅ Configured: {info['configured']}")
    
    # Show available biotypes
    data_avail = info.get("data_available", {})
    if "biotypes" in data_avail:
        biotypes = data_avail["biotypes"]
        print(f"✅ Gene biotypes available: {len(biotypes)}")
        print(f"   Types: {', '.join(biotypes)}")
    
    # Show gene counts
    if "gene_counts" in data_avail:
        gene_counts = data_avail["gene_counts"]
        total_genes = sum(gene_counts.values())
        print(f"✅ Total genes: {total_genes:,}")
        
        print("   Breakdown:")
        for biotype, count in sorted(gene_counts.items()):
            percentage = (count / total_genes) * 100
            print(f"     {biotype}: {count:,} genes ({percentage:.1f}%)")
    
    # Show chromosome data
    if "chromosomes" in data_avail:
        chroms = data_avail["chromosomes"]
        print(f"✅ Chromosomes: {len(chroms)} available")
        print(f"   Available: {', '.join(sorted(chroms)[:10])}{'...' if len(chroms) > 10 else ''}")
    
    print()

def gene_exploration_examples():
    """Demonstrate gene exploration utilities"""
    print("=== Gene Exploration Examples ===")
    
    available = list_available_organisms()
    if not available:
        print("No organisms configured. Run setup first.")
        return
    
    organism = available[0]
    print(f"Exploring genes in: {organism}")
    
    # List gene biotypes
    biotypes = list_gene_biotypes(organism)
    print(f"✅ Available biotypes: {biotypes}")
    
    # Count genes for each biotype
    counts = count_genes(organism)
    print(f"✅ Gene counts: {counts}")
    
    # Focus on protein-coding genes
    if 'protein_coding' in biotypes:
        print(f"\n--- Protein-coding genes ---")
        pc_count = counts.get('protein_coding', 0)
        print(f"Total protein-coding genes: {pc_count:,}")
        
        # Get a sample of gene names
        gene_sample = get_gene_list(organism, 'protein_coding', limit=10)
        print(f"Sample genes: {gene_sample}")
        
        # Search for specific genes
        print(f"\n--- Gene search examples ---")
        
        # Search for genes containing 'TP'
        tp_genes = search_genes(organism, 'TP', biotype='protein_coding', limit=5)
        print(f"Genes containing 'TP': {[g['gene_name'] for g in tp_genes]}")
        
        # Search for genes starting with 'K'
        k_genes = search_genes(organism, 'K', biotype='protein_coding', limit=5)
        print(f"Genes starting with 'K': {[g['gene_name'] for g in k_genes]}")
    
    print()

def comprehensive_summary_example():
    """Show comprehensive data summary"""
    print("=== Comprehensive Data Summary ===")
    
    # Get programmatic summary
    summary = data_summary()
    
    print("Programmatic summary:")
    print(f"  Supported organisms: {summary['supported_organisms']}")
    print(f"  Configured organisms: {summary['configured_organisms']}")
    print(f"  Totals: {summary['totals']}")
    
    print("\nFormatted summary:")
    print("-" * 50)
    
    # Print formatted summary
    print_data_summary()

def search_examples():
    """Demonstrate various search patterns"""
    print("=== Gene Search Examples ===")
    
    available = list_available_organisms()
    if not available:
        print("No organisms configured. Run setup first.")
        return
    
    organism = available[0]
    
    search_patterns = [
        ("KRAS", "Exact gene name"),
        ("K", "Genes starting with K"),
        ("ase", "Genes containing 'ase' (enzymes)"),
        ("P53", "Genes containing 'P53'"),
        ("HOX", "HOX genes (development)")
    ]
    
    for pattern, description in search_patterns:
        print(f"\n--- {description} ---")
        results = search_genes(organism, pattern, limit=5)
        
        if results:
            for i, result in enumerate(results, 1):
                print(f"  {i}. {result['gene_name']} ({result['biotype']})")
        else:
            print(f"  No genes found matching '{pattern}'")

def data_exploration_workflow():
    """Complete workflow for exploring new genomics data"""
    print("=== Complete Data Exploration Workflow ===")
    
    print("Step 1: Check organism support")
    supported = list_supported_organisms()
    available = list_available_organisms()
    
    print(f"  Supported: {supported}")
    print(f"  Configured: {available}")
    
    if not available:
        print("\n  🔄 To set up data:")
        print("    seqmat setup --path /your/data/path --organism hg38")
        return
    
    for organism in available:
        print(f"\nStep 2: Explore {organism} data")
        
        # Get overview
        info = get_organism_info(organism)
        if "error" in info:
            continue
        
        biotypes = list_gene_biotypes(organism)
        counts = count_genes(organism)
        
        print(f"  Biotypes: {len(biotypes)}")
        print(f"  Total genes: {sum(counts.values()):,}")
        
        # Find most common biotype
        if counts:
            most_common = max(counts.items(), key=lambda x: x[1])
            print(f"  Most common: {most_common[0]} ({most_common[1]:,} genes)")
            
            # Sample genes from most common biotype
            sample_genes = get_gene_list(organism, most_common[0], limit=3)
            print(f"  Sample genes: {sample_genes}")
    
    print(f"\nStep 3: Summary")
    summary = data_summary()
    totals = summary['totals']
    print(f"  Ready to analyze {totals['genes']:,} genes across {totals['organisms']} organisms!")

def main():
    """Run all examples"""
    print("🔍 SeqMat Data Inspection Examples\n")
    
    examples = [
        basic_organism_info,
        explore_organism_data,
        gene_exploration_examples,
        search_examples,
        comprehensive_summary_example,
        data_exploration_workflow
    ]
    
    for example_func in examples:
        try:
            example_func()
        except Exception as e:
            print(f"⚠️  Example {example_func.__name__} encountered an error: {e}")
            print("This is likely because no genomics data has been set up yet.")
            print("Run: seqmat setup --path /your/data/path --organism hg38")
            print()
            continue
    
    print("=" * 60)
    print("🎯 Data Inspection Features Summary:")
    print("✅ Organism support checking")
    print("✅ Data availability inspection") 
    print("✅ Gene biotype enumeration")
    print("✅ Gene counting and listing")
    print("✅ Flexible gene search")
    print("✅ Comprehensive data summaries")
    print("✅ Both Python API and CLI access")
    print("\n🚀 Ready for genomics data exploration!")

if __name__ == "__main__":
    main()