"""
Basic SeqMat Usage Examples

This script demonstrates the core functionality of the SeqMat library.
"""

from seqmat import SeqMat
import numpy as np

def basic_sequence_operations():
    """Demonstrate basic sequence creation and manipulation"""
    print("=== Basic Sequence Operations ===")
    
    # Create a simple sequence
    seq = SeqMat("ATCGATCGATCG", name="example_sequence")
    print(f"Original sequence: {seq.seq}")
    print(f"Length: {len(seq)} bp")
    print(f"Name: {seq.name}")
    print()
    
    # Create sequence with genomic coordinates
    indices = np.arange(1000, 1012)
    genomic_seq = SeqMat("ATCGATCGATCG", indices=indices, name="chr1:1000-1011")
    print(f"Genomic sequence: {genomic_seq.seq}")
    print(f"Coordinates: {genomic_seq.index[0]}-{genomic_seq.index[-1]}")
    print()

def mutation_examples():
    """Demonstrate mutation operations"""
    print("=== Mutation Examples ===")
    
    seq = SeqMat("ATCGATCGATCG", name="mutation_demo")
    print(f"Original: {seq.seq}")
    
    # Single nucleotide polymorphism
    seq.apply_mutations((5, "T", "A"))
    print(f"After SNP (T->A at pos 5): {seq.seq}")
    
    # Insertion
    seq.apply_mutations((8, "-", "GGG"))
    print(f"After insertion (GGG at pos 8): {seq.seq}")
    
    # Deletion
    seq.apply_mutations((3, "GA", "-"))
    print(f"After deletion (GA at pos 3): {seq.seq}")
    
    print(f"Total mutations applied: {len(seq.mutations)}")
    print(f"Modified positions: {sorted(seq.mutated_positions)}")
    print()

def sequence_transformations():
    """Demonstrate sequence transformations"""
    print("=== Sequence Transformations ===")
    
    seq = SeqMat("ATCGATCG", name="transform_demo")
    print(f"Original:     {seq.seq}")
    
    # Complement
    comp = seq.complement()
    print(f"Complement:   {comp.seq}")
    
    # Reverse complement
    rev_comp = seq.clone()
    rev_comp.reverse_complement()
    print(f"Rev comp:     {rev_comp.seq}")
    print(f"Reversed:     {rev_comp.rev}")
    print()

def slicing_and_indexing():
    """Demonstrate slicing and indexing operations"""
    print("=== Slicing and Indexing ===")
    
    # Create sequence with specific coordinates
    indices = np.arange(100, 112)
    seq = SeqMat("ATCGATCGATCG", indices=indices, name="slice_demo")
    print(f"Full sequence: {seq.seq} (positions {indices[0]}-{indices[-1]})")
    
    # Slice by position range
    subseq = seq[103:108]
    print(f"Slice [103:108]: {subseq.seq}")
    print(f"Slice coordinates: {subseq.index}")
    
    # Access single position
    base_info = seq[105]
    print(f"Base at position 105: {base_info['nt'].decode()}")
    print()

def region_removal():
    """Demonstrate region removal (e.g., splicing)"""
    print("=== Region Removal (Splicing) ===")
    
    seq = SeqMat("ATCGATCGATCGATCG", indices=np.arange(1, 17), name="splicing_demo")
    print(f"Pre-mRNA:  {seq.seq}")
    
    # Define introns to remove
    introns = [(4, 7), (11, 13)]
    print(f"Introns to remove: {introns}")
    
    # Splice out introns
    mature = seq.remove_regions(introns)
    print(f"Mature mRNA: {mature.seq}")
    print(f"Length change: {len(seq)} -> {len(mature)} bp")
    print()

def advanced_mutations():
    """Demonstrate advanced mutation features"""
    print("=== Advanced Mutations ===")
    
    seq = SeqMat("ATCGATCGATCGATCG", name="advanced_demo")
    print(f"Original: {seq.seq}")
    
    # Multiple mutations at once
    mutations = [
        (2, "T", "A"),      # SNP
        (6, "-", "AAA"),    # Insertion
        (10, "CGA", "-"),   # Deletion
        (14, "TC", "GG")    # Complex substitution
    ]
    
    seq.apply_mutations(mutations)
    print(f"After multiple mutations: {seq.seq}")
    
    # Show mutation history
    print("\nMutation history:")
    for i, mut in enumerate(seq.mutations, 1):
        print(f"  {i}. {mut['type'].upper()} at pos {mut['pos']}: "
              f"{mut['ref']} -> {mut['alt']}")
    
    # Try conflicting mutations (will be detected and skipped)
    print("\nTrying conflicting mutations...")
    conflicting = [
        (5, "A", "T"),
        (5, "T", "G")  # Conflicts with previous
    ]
    seq.apply_mutations(conflicting)
    print()

def export_operations():
    """Demonstrate export functionality"""
    print("=== Export Operations ===")
    
    seq = SeqMat("ATCGATCGATCG", name="export_demo")
    seq.apply_mutations([(5, "T", "A"), (8, "-", "GG")])
    
    # Export as FASTA
    fasta = seq.to_fasta(wrap=40)
    print("FASTA format:")
    print(fasta)
    print()
    
    # Export as dictionary
    data = seq.to_dict()
    print("Dictionary format:")
    for key, value in data.items():
        print(f"  {key}: {value}")
    print()

def main():
    """Run all examples"""
    basic_sequence_operations()
    mutation_examples()
    sequence_transformations()
    slicing_and_indexing()
    region_removal()
    advanced_mutations()
    export_operations()
    
    print("=== Summary ===")
    print("This script demonstrated:")
    print("- Basic sequence creation and manipulation")
    print("- Mutation operations (SNPs, insertions, deletions)")
    print("- Sequence transformations (complement, reverse)")
    print("- Slicing and indexing by genomic coordinates")
    print("- Region removal for splicing analysis")
    print("- Advanced mutation features and conflict detection")
    print("- Export operations (FASTA, dictionary)")

if __name__ == "__main__":
    main()