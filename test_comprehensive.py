#!/usr/bin/env python3
"""
Comprehensive test script for all SeqMat, Gene, and Transcript functionalities
"""

import seqmat
from seqmat import SeqMat, Gene, Transcript
from seqmat import (
    setup_genomics_data, list_available_organisms, list_supported_organisms,
    get_organism_info, list_gene_biotypes, count_genes, get_gene_list,
    data_summary, search_genes, get_all_genes
)
import time

def print_section(title):
    """Print a formatted section header"""
    print(f"\n{'='*60}")
    print(f"{title.center(60)}")
    print(f"{'='*60}\n")

def test_config_utilities():
    """Test configuration and utility functions"""
    print_section("Testing Configuration and Utility Functions")
    
    # Test 1: List supported organisms
    print("1. Supported organisms:")
    supported = list_supported_organisms()
    for org in supported:
        print(f"   - {org}")
    
    # Test 2: List configured organisms
    print("\n2. Configured organisms:")
    configured = list_available_organisms()
    for org in configured:
        print(f"   - {org}")
    
    # Test 3: Get default organism
    default_org = seqmat.get_default_organism()
    print(f"\n3. Default organism: {default_org}")
    
    # Test 4: Get organism info
    print(f"\n4. Organism info for {default_org}:")
    info = get_organism_info(default_org)
    if 'data_available' in info:
        print(f"   Biotypes: {info['data_available'].get('biotypes', [])}")
        print(f"   Gene counts: {info['data_available'].get('gene_counts', {})}")
        print(f"   Chromosomes: {len(info['data_available'].get('chromosomes', []))}")
    
    # Test 5: List gene biotypes
    print(f"\n5. Gene biotypes for {default_org}:")
    biotypes = list_gene_biotypes(default_org)
    for bt in biotypes:
        print(f"   - {bt}")
    
    # Test 6: Count genes
    print(f"\n6. Gene counts by biotype:")
    gene_counts = count_genes(default_org)
    for bt, count in gene_counts.items():
        print(f"   {bt}: {count}")
    
    # Test 7: Get gene list (limited)
    print(f"\n7. First 5 protein-coding genes:")
    genes = get_gene_list(default_org, 'protein_coding', limit=5)
    for gene in genes:
        print(f"   - {gene}")
    
    return default_org

def test_gene_search_functions(organism):
    """Test gene search and retrieval functions"""
    print_section("Testing Gene Search Functions")
    
    # Test 1: Search genes
    print("1. Searching for genes containing 'KRAS':")
    results = search_genes(organism, 'KRAS', limit=5)
    for result in results:
        print(f"   - {result['gene_name']} ({result['gene_id']}) - {result['biotype']}")
    
    # Test 2: Get all genes (limited)
    print("\n2. Getting all genes (first 10):")
    all_genes = get_all_genes(organism)
    print(f"   Total genes: {len(all_genes)}")
    for gene in all_genes[:10]:
        print(f"   - {gene['gene_name']} ({result['gene_id']}) - {gene['biotype']}")
    
    # Test 3: Get genes by biotype
    print("\n3. Getting lncRNA genes (first 5):")
    lncrna_genes = get_all_genes(organism, biotype='lncRNA')
    print(f"   Total lncRNA genes: {len(lncrna_genes)}")
    for gene in lncrna_genes[:5]:
        print(f"   - {gene['gene_name']} ({gene['gene_id']})")
    
    return results[0]['gene_name'] if results else None

def test_seqmat_basic():
    """Test basic SeqMat functionality"""
    print_section("Testing Basic SeqMat Functionality")
    
    # Test 1: Create SeqMat from string
    print("1. Creating SeqMat from string:")
    seq = SeqMat("ATCGATCGATCG")
    print(f"   Sequence: {seq.seq}")
    print(f"   Length: {len(seq)}")
    print(f"   String representation: {str(seq)}")
    
    # Test 2: Complement
    print("\n2. Complement:")
    comp = seq.complement()
    print(f"   Original: {seq.seq}")
    print(f"   Complement: {comp.seq}")
    
    # Test 3: Reverse complement
    print("\n3. Reverse complement:")
    seq2 = SeqMat("ATCGATCGATCG")
    seq2.reverse_complement()
    print(f"   Original: ATCGATCGATCG")
    print(f"   Rev comp: {seq2.seq}")
    
    # Test 4: Clone
    print("\n4. Cloning:")
    clone = seq.clone()
    clone.seq = "AAAA"
    print(f"   Original: {seq.seq}")
    print(f"   Clone (modified): {clone.seq}")
    
    # Test 5: Subsequence
    print("\n5. Subsequence extraction:")
    subseq = seq[3:9]
    print(f"   Full sequence: {seq.seq}")
    print(f"   Subsequence [3:9]: {subseq.seq}")
    
    # Test 6: SeqMat with metadata
    print("\n6. SeqMat with metadata:")
    seq_meta = SeqMat("ATCG", name="test_seq", organism="hg38")
    print(f"   Sequence: {seq_meta.seq}")
    print(f"   Name: {seq_meta.name}")
    print(f"   Organism: {seq_meta.organism}")

def test_seqmat_mutations():
    """Test SeqMat mutation functionality"""
    print_section("Testing SeqMat Mutations")
    
    # Test 1: Single SNP
    print("1. Single SNP mutation:")
    seq = SeqMat("ATCGATCGATCG")
    mutations = [(3, "G", "C")]
    seq.apply_mutations(mutations)
    print(f"   Original: ATCGATCGATCG")
    print(f"   Mutated:  {seq.seq}")
    print(f"   Mutations: {seq.mutations}")
    
    # Test 2: Multiple mutations
    print("\n2. Multiple mutations (SNP, insertion, deletion):")
    seq = SeqMat("ATCGATCGATCGATCG")
    mutations = [
        (3, "G", "C"),      # SNP
        (6, "-", "AAA"),    # Insertion
        (12, "ATC", "-")    # Deletion
    ]
    seq.apply_mutations(mutations)
    print(f"   Original: ATCGATCGATCGATCG")
    print(f"   Mutated:  {seq.seq}")
    print(f"   Number of mutations: {len(seq.mutations)}")
    
    # Test 3: Mutation with quiet mode
    print("\n3. Mutations in quiet mode:")
    seq = SeqMat("ATCGATCGATCG")
    seq.apply_mutations([(3, "G", "T")], quiet=True)
    print(f"   Mutated sequence: {seq.seq}")
    print(f"   No position tracking in quiet mode")
    
    # Test 4: Chained mutations
    print("\n4. Chained mutations:")
    seq = SeqMat("ATCGATCGATCG")
    seq.apply_mutations([(3, "G", "C")])
    print(f"   After first mutation: {seq.seq}")
    seq.apply_mutations([(6, "C", "T")])
    print(f"   After second mutation: {seq.seq}")
    print(f"   Total mutations: {len(seq.mutations)}")
    
    # Test 5: Strand-aware mutations
    print("\n5. Strand-aware mutations:")
    seq = SeqMat("ATCGATCGATCG", rev=True)
    mutations = [(3, "G", "C")]
    seq.apply_mutations(mutations)
    print(f"   Reverse strand sequence: {seq.seq}")
    print(f"   Mutations applied considering reverse strand")

def test_seqmat_fasta(organism):
    """Test SeqMat FASTA loading functionality"""
    print_section("Testing SeqMat FASTA Loading")
    
    # Get available chromosomes
    info = get_organism_info(organism)
    chromosomes = info.get('data_available', {}).get('chromosomes', [])
    
    if not chromosomes:
        print("No chromosomes available for testing")
        return
    
    # Use first available chromosome
    chrom = chromosomes[0]
    print(f"Using chromosome: {chrom}")
    
    # Test 1: Load sequence region
    print("\n1. Loading sequence region:")
    seq = SeqMat.from_fasta(
        genome=organism,
        chrom=chrom,
        start=1000000,
        end=1000100
    )
    print(f"   Loaded {len(seq)}bp from {chrom}:1000000-1000100")
    print(f"   First 50bp: {seq.seq[:50]}...")
    
    # Test 2: Load with mutations
    print("\n2. Loading with mutations:")
    mutations = [(1000010, "N", "A"), (1000020, "N", "T")]
    seq_mut = SeqMat.from_fasta(
        genome=organism,
        chrom=chrom,
        start=1000000,
        end=1000100,
        mutations=mutations
    )
    print(f"   Loaded with {len(mutations)} mutations")
    print(f"   Applied mutations: {len(seq_mut.mutations)}")
    
    # Test 3: Load reverse strand
    print("\n3. Loading reverse strand:")
    seq_rev = SeqMat.from_fasta(
        genome=organism,
        chrom=chrom,
        start=1000000,
        end=1000100,
        rev=True
    )
    print(f"   Reverse strand loaded")
    print(f"   First 50bp: {seq_rev.seq[:50]}...")

def test_gene_functionality(organism):
    """Test Gene class functionality"""
    print_section("Testing Gene Class")
    
    # Find a test gene
    genes = get_gene_list(organism, 'protein_coding', limit=10)
    if not genes:
        print("No genes available for testing")
        return None
    
    test_gene_name = genes[0]
    print(f"Using test gene: {test_gene_name}")
    
    # Test 1: Load gene from file
    print("\n1. Loading gene from file:")
    gene = Gene.from_file(test_gene_name, organism=organism)
    if not gene:
        print("   Failed to load gene")
        return None
    
    print(f"   Gene name: {gene.gene_name}")
    print(f"   Gene ID: {gene.gene_id}")
    print(f"   Chromosome: {gene.chrm}")
    print(f"   Reverse strand: {gene.rev}")
    print(f"   Number of transcripts: {len(gene.transcripts)}")
    
    # Test 2: String representations
    print("\n2. String representations:")
    print(f"   repr: {repr(gene)}")
    print(f"   str: {str(gene)}")
    
    # Test 3: Iterate over transcripts
    print("\n3. Iterating over transcripts (first 3):")
    for i, transcript in enumerate(gene):
        if i >= 3:
            break
        print(f"   - {transcript.transcript_id} ({transcript.transcript_biotype})")
    
    # Test 4: Get specific transcript
    print("\n4. Getting specific transcript:")
    tid = list(gene.transcripts.keys())[0] if gene.transcripts else None
    if tid:
        transcript = gene[tid]
        print(f"   Retrieved transcript: {transcript.transcript_id}")
    
    # Test 5: Get primary transcript
    print("\n5. Getting primary transcript:")
    primary = gene.transcript()
    if primary:
        print(f"   Primary transcript: {primary.transcript_id}")
        print(f"   Biotype: {primary.transcript_biotype}")
    
    # Test 6: Splice sites
    print("\n6. Analyzing splice sites:")
    acceptors, donors = gene.splice_sites()
    print(f"   Unique acceptor sites: {len(acceptors)}")
    print(f"   Unique donor sites: {len(donors)}")
    
    # Test 7: Test non-protein-coding gene
    print("\n7. Testing non-protein-coding gene:")
    lncrna_genes = get_all_genes(organism, biotype='lncRNA')
    if lncrna_genes:
        lncrna_name = lncrna_genes[0]['gene_name']
        lncrna_gene = Gene.from_file(lncrna_name, organism=organism)
        if lncrna_gene:
            print(f"   Loaded lncRNA gene: {lncrna_gene.gene_name}")
            print(f"   Transcripts: {len(lncrna_gene.transcripts)}")
    
    return gene

def test_transcript_functionality(gene):
    """Test Transcript class functionality"""
    print_section("Testing Transcript Class")
    
    if not gene or not gene.transcripts:
        print("No gene available for transcript testing")
        return
    
    # Get a transcript
    transcript = gene.transcript()
    if not transcript:
        print("No primary transcript available")
        return
    
    print(f"Testing transcript: {transcript.transcript_id}")
    
    # Test 1: Basic properties
    print("\n1. Basic transcript properties:")
    print(f"   Transcript ID: {transcript.transcript_id}")
    print(f"   Biotype: {transcript.transcript_biotype}")
    print(f"   Start: {transcript.transcript_start}")
    print(f"   End: {transcript.transcript_end}")
    print(f"   Length: {len(transcript)}")
    print(f"   Reverse strand: {transcript.rev}")
    print(f"   Chromosome: {transcript.chrm}")
    
    # Test 2: Exon information
    print("\n2. Exon information:")
    print(f"   Number of exons: {len(transcript.exons)}")
    if transcript.exons:
        print(f"   First exon: {transcript.exons[0]}")
        print(f"   Last exon: {transcript.exons[-1]}")
    
    # Test 3: Protein coding features
    print("\n3. Protein coding features:")
    print(f"   Is protein coding: {transcript.protein_coding}")
    if transcript.protein_coding:
        print(f"   Protein ID: {transcript.protein_id}")
        print(f"   TIS (Translation Initiation Site): {transcript.TIS}")
        print(f"   TTS (Translation Termination Site): {transcript.TTS}")
    
    # Test 4: Splice sites
    print("\n4. Splice sites:")
    if hasattr(transcript, 'acceptors') and transcript.acceptors:
        print(f"   Acceptor sites: {len(transcript.acceptors)}")
        print(f"   First acceptor: {transcript.acceptors[0]}")
    if hasattr(transcript, 'donors') and transcript.donors:
        print(f"   Donor sites: {len(transcript.donors)}")
        print(f"   First donor: {transcript.donors[0]}")
    
    # Test 5: Generate mature mRNA
    print("\n5. Generating mature mRNA:")
    try:
        transcript.generate_mature_mrna()
        if transcript.mature_mrna:
            print(f"   Mature mRNA length: {len(transcript.mature_mrna)}bp")
            print(f"   First 50bp: {transcript.mature_mrna[:50]}...")
        else:
            print("   Failed to generate mature mRNA")
    except Exception as e:
        print(f"   Error generating mature mRNA: {e}")
    
    # Test 6: CDS sequence
    print("\n6. CDS (Coding Sequence):")
    if transcript.protein_coding:
        try:
            cds = transcript.cds()
            if cds:
                print(f"   CDS length: {len(cds)}bp")
                print(f"   First 30bp: {cds[:30]}...")
                # Check if divisible by 3
                print(f"   Valid CDS length (divisible by 3): {len(cds) % 3 == 0}")
        except Exception as e:
            print(f"   Error getting CDS: {e}")
    
    # Test 7: Conservation data
    print("\n7. Conservation data:")
    print(f"   Conservation available: {transcript.cons_available}")
    
    # Test 8: Primary transcript flag
    print("\n8. Primary transcript:")
    print(f"   Is primary transcript: {transcript.primary_transcript}")
    
    # Test 9: String representations
    print("\n9. String representations:")
    print(f"   repr: {repr(transcript)}")
    print(f"   str: {str(transcript)}")

def test_data_summary():
    """Test data summary functionality"""
    print_section("Testing Data Summary")
    
    # Get comprehensive summary
    summary = data_summary()
    
    print("1. Summary statistics:")
    print(f"   Supported organisms: {summary['totals']['organisms']}")
    print(f"   Total biotypes: {summary['totals']['biotypes']}")
    print(f"   Total genes: {summary['totals']['genes']}")
    
    print("\n2. Configured organisms:")
    for org in summary['configured_organisms']:
        print(f"   - {org}")

def run_performance_test(organism):
    """Run basic performance tests"""
    print_section("Performance Tests")
    
    # Test 1: Gene loading speed
    print("1. Gene loading performance:")
    genes = get_gene_list(organism, 'protein_coding', limit=10)
    if genes:
        start_time = time.time()
        for gene_name in genes[:5]:
            gene = Gene.from_file(gene_name, organism=organism)
        end_time = time.time()
        avg_time = (end_time - start_time) / 5
        print(f"   Average time per gene load: {avg_time:.3f}s")
    
    # Test 2: Sequence manipulation speed
    print("\n2. Sequence manipulation performance:")
    seq = SeqMat("A" * 10000)  # 10kb sequence
    
    start_time = time.time()
    seq.reverse_complement()
    rc_time = time.time() - start_time
    print(f"   Reverse complement of 10kb: {rc_time:.3f}s")
    
    start_time = time.time()
    mutations = [(i*100, "A", "T") for i in range(10)]
    seq.apply_mutations(mutations)
    mut_time = time.time() - start_time
    print(f"   Apply 10 mutations to 10kb: {mut_time:.3f}s")
    
    # Test 3: Search performance
    print("\n3. Search performance:")
    start_time = time.time()
    results = search_genes(organism, "A", limit=100)
    search_time = time.time() - start_time
    print(f"   Search for 'A' (100 results): {search_time:.3f}s")
    print(f"   Results found: {len(results)}")

def main():
    """Main test function"""
    print("=" * 60)
    print("COMPREHENSIVE SEQMAT TEST SUITE".center(60))
    print("=" * 60)
    
    # Test configuration utilities
    organism = test_config_utilities()
    
    # Test gene search functions
    test_gene_name = test_gene_search_functions(organism)
    
    # Test SeqMat basic functionality
    test_seqmat_basic()
    
    # Test SeqMat mutations
    test_seqmat_mutations()
    
    # Test SeqMat FASTA loading
    test_seqmat_fasta(organism)
    
    # Test Gene functionality
    gene = test_gene_functionality(organism)
    
    # Test Transcript functionality
    if gene:
        test_transcript_functionality(gene)
    
    # Test data summary
    test_data_summary()
    
    # Run performance tests
    run_performance_test(organism)
    
    print("\n" + "=" * 60)
    print("ALL TESTS COMPLETED".center(60))
    print("=" * 60)

if __name__ == "__main__":
    main()