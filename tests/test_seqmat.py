"""Tests for SeqMat class"""
import pytest
import numpy as np
from seqmat import SeqMat


class TestSeqMat:
    """Test cases for SeqMat functionality"""
    
    def test_creation_from_string(self):
        """Test creating SeqMat from string"""
        seq = SeqMat("ATCGATCG", name="test_seq")
        assert len(seq) == 8
        assert seq.seq == "ATCGATCG"
        assert seq.name == "test_seq"
    
    def test_creation_with_indices(self):
        """Test creating SeqMat with custom indices"""
        indices = np.array([100, 101, 102, 103])
        seq = SeqMat("ATCG", indices=indices)
        assert np.array_equal(seq.index, indices)
    
    def test_snp_mutation(self):
        """Test single nucleotide polymorphism"""
        seq = SeqMat("ATCGATCG")
        seq.apply_mutations((3, "C", "G"))
        assert seq.seq == "ATGGATCG"
        assert 3 in seq.mutated_positions
    
    def test_insertion_mutation(self):
        """Test insertion mutation"""
        seq = SeqMat("ATCGATCG")
        seq.apply_mutations((4, "-", "AAA"))
        # Length should increase
        assert len(seq) == 11
        assert "AAA" in seq.seq
    
    def test_deletion_mutation(self):
        """Test deletion mutation"""
        seq = SeqMat("ATCGATCG")
        seq.apply_mutations((3, "CGA", "-"))
        assert len(seq) == 5
        # Check that deleted bases are not in valid sequence
        assert seq.seq == "ATTCG"
    
    def test_multiple_mutations(self):
        """Test applying multiple mutations"""
        seq = SeqMat("ATCGATCGATCG")
        mutations = [
            (2, "T", "A"),      # SNP
            (6, "-", "GGG"),    # Insertion
            (10, "AT", "-")     # Deletion
        ]
        seq.apply_mutations(mutations)
        assert len(seq.mutations) == 3
        assert len(seq.mutated_positions) >= 3
    
    def test_complement(self):
        """Test complement operation"""
        seq = SeqMat("ATCG")
        comp = seq.complement()
        assert comp.seq == "TAGC"
    
    def test_reverse_complement(self):
        """Test reverse complement operation"""
        seq = SeqMat("ATCG")
        seq.reverse_complement()
        assert seq.seq == "CGAT"
        assert seq.rev == True
    
    def test_slicing(self):
        """Test sequence slicing"""
        seq = SeqMat("ATCGATCGATCG", indices=np.arange(100, 112))
        subseq = seq[103:107]
        assert len(subseq) == 4
        assert subseq.index[0] == 103
        assert subseq.index[-1] == 106
    
    def test_remove_regions(self):
        """Test removing regions from sequence"""
        seq = SeqMat("ATCGATCGATCG", indices=np.arange(1, 13))
        regions = [(3, 5), (8, 9)]
        new_seq = seq.remove_regions(regions)
        # Should have removed positions 3-5 and 8-9
        assert len(new_seq) < len(seq)
    
    def test_mutation_validation(self):
        """Test mutation conflict detection"""
        seq = SeqMat("ATCGATCG")
        # Overlapping mutations should be caught
        mutations = [
            (3, "CGA", "TTT"),
            (4, "GAT", "AAA")  # Overlaps with previous
        ]
        result = seq._validate_mutation_batch(mutations)
        assert result == False
    
    def test_clone(self):
        """Test cloning functionality"""
        seq = SeqMat("ATCG")
        seq.apply_mutations((2, "T", "A"))
        
        clone = seq.clone()
        clone.apply_mutations((4, "G", "C"))
        
        # Original should be unchanged
        assert seq.seq == "AACG"
        assert clone.seq == "AACC"
        assert len(seq.mutations) == 1
        assert len(clone.mutations) == 2
    
    def test_to_fasta(self):
        """Test FASTA export"""
        seq = SeqMat("ATCGATCG", name="test_sequence")
        fasta = seq.to_fasta()
        assert fasta.startswith(">test_sequence")
        assert "ATCGATCG" in fasta
    
    def test_to_dict(self):
        """Test dictionary conversion"""
        seq = SeqMat("ATCG", name="test")
        seq.apply_mutations((2, "T", "A"))
        
        data = seq.to_dict()
        assert data['name'] == 'test'
        assert data['sequence'] == 'AACG'
        assert data['length'] == 4
        assert len(data['mutations']) == 1