"""Tests for Gene class"""
import pytest
from unittest.mock import patch, MagicMock
from seqmat import Gene


class TestGene:
    """Test cases for Gene functionality"""
    
    def test_gene_initialization(self):
        """Test Gene object creation"""
        gene = Gene(
            gene_name="KRAS",
            gene_id="ENSG00000133703", 
            rev=False,
            chrm="12",
            transcripts={"ENST00000311936": {"transcript_id": "ENST00000311936"}},
            organism="hg38"
        )
        
        assert gene.gene_name == "KRAS"
        assert gene.gene_id == "ENSG00000133703"
        assert gene.chrm == "12"
        assert len(gene) == 1
    
    def test_gene_string_representations(self):
        """Test __str__ and __repr__ methods"""
        gene = Gene("KRAS", "ENSG00000133703", False, "12", {})
        
        assert "KRAS" in str(gene)
        assert "ENSG00000133703" in str(gene)
        assert repr(gene) == "Gene(KRAS)"
    
    def test_gene_iteration(self):
        """Test iterating over transcripts"""
        transcripts = {
            "ENST1": {
                "transcript_id": "ENST1",
                "transcript_start": 100,
                "transcript_end": 200,
                "rev": False,
                "chrm": "12"
            },
            "ENST2": {
                "transcript_id": "ENST2", 
                "transcript_start": 150,
                "transcript_end": 250,
                "rev": False,
                "chrm": "12"
            }
        }
        
        gene = Gene("TEST", "ENSG123", False, "12", transcripts)
        transcript_ids = [t.transcript_id for t in gene]
        
        assert "ENST1" in transcript_ids
        assert "ENST2" in transcript_ids
        assert len(transcript_ids) == 2
    
    def test_gene_getitem(self):
        """Test accessing transcripts by ID"""
        transcript_data = {
            "transcript_id": "ENST1",
            "transcript_start": 100,
            "transcript_end": 200,
            "rev": False,
            "chrm": "12"
        }
        
        gene = Gene("TEST", "ENSG123", False, "12", {"ENST1": transcript_data})
        
        # Valid transcript
        transcript = gene["ENST1"]
        assert transcript is not None
        assert transcript.transcript_id == "ENST1"
        
        # Invalid transcript
        invalid = gene["NONEXISTENT"]
        assert invalid is None
    
    def test_splice_sites(self):
        """Test splice site aggregation"""
        transcripts = {
            "ENST1": {
                "transcript_id": "ENST1",
                "acceptors": [100, 200],
                "donors": [150, 250],
                "transcript_start": 50,
                "transcript_end": 300,
                "rev": False,
                "chrm": "12"
            },
            "ENST2": {
                "transcript_id": "ENST2",
                "acceptors": [100, 220],  # 100 is shared
                "donors": [180, 250],     # 250 is shared
                "transcript_start": 50,
                "transcript_end": 300,
                "rev": False,
                "chrm": "12"
            }
        }
        
        gene = Gene("TEST", "ENSG123", False, "12", transcripts)
        acceptors, donors = gene.splice_sites()
        
        # Should count shared sites
        assert acceptors[100] == 2  # Shared acceptor
        assert donors[250] == 2     # Shared donor
        assert len(acceptors) == 3  # 100, 200, 220
        assert len(donors) == 3     # 150, 180, 250
    
    def test_primary_transcript_selection(self):
        """Test primary transcript identification"""
        transcripts = {
            "ENST1": {
                "transcript_id": "ENST1",
                "primary_transcript": False,
                "transcript_biotype": "protein_coding",
                "transcript_start": 100,
                "transcript_end": 200,
                "rev": False,
                "chrm": "12"
            },
            "ENST2": {
                "transcript_id": "ENST2",
                "primary_transcript": True, 
                "transcript_biotype": "protein_coding",
                "transcript_start": 150,
                "transcript_end": 250,
                "rev": False,
                "chrm": "12"
            }
        }
        
        gene = Gene("TEST", "ENSG123", False, "12", transcripts)
        primary_id = gene.primary_transcript
        
        assert primary_id == "ENST2"
        
        # Test getting primary transcript object
        primary = gene.transcript()
        assert primary.transcript_id == "ENST2"
    
    def test_primary_transcript_fallback(self):
        """Test fallback to protein-coding when no primary marked"""
        transcripts = {
            "ENST1": {
                "transcript_id": "ENST1",
                "transcript_biotype": "nonsense_mediated_decay",
                "transcript_start": 100,
                "transcript_end": 200,
                "rev": False,
                "chrm": "12"
            },
            "ENST2": {
                "transcript_id": "ENST2",
                "transcript_biotype": "protein_coding",
                "transcript_start": 150,
                "transcript_end": 250,
                "rev": False,
                "chrm": "12"
            }
        }
        
        gene = Gene("TEST", "ENSG123", False, "12", transcripts)
        primary_id = gene.primary_transcript
        
        assert primary_id == "ENST2"  # Should pick protein-coding
    
    @patch('seqmat.gene.get_organism_config')
    @patch('seqmat.gene.unload_pickle')
    def test_from_file_success(self, mock_unload, mock_config):
        """Test successful gene loading from file"""
        # Mock configuration
        mock_config.return_value = {
            'MRNA_PATH': '/mock/path'
        }
        
        # Mock file data
        gene_data = {
            'gene_name': 'KRAS',
            'gene_id': 'ENSG00000133703',
            'rev': False,
            'chrm': '12',
            'transcripts': {}
        }
        mock_unload.return_value = gene_data
        
        # Mock file existence
        with patch('pathlib.Path.glob') as mock_glob:
            mock_glob.return_value = ['mock_file.pkl']
            
            gene = Gene.from_file('KRAS', 'hg38')
            
            assert gene is not None
            assert gene.gene_name == 'KRAS'
            assert gene.gene_id == 'ENSG00000133703'
    
    @patch('seqmat.gene.get_organism_config')
    def test_from_file_not_configured(self, mock_config):
        """Test gene loading when organism not configured"""
        mock_config.side_effect = ValueError("Organism not configured")
        
        gene = Gene.from_file('KRAS', 'hg38')
        assert gene is None
    
    @patch('seqmat.gene.get_organism_config')
    def test_from_file_not_found(self, mock_config):
        """Test gene loading when file not found"""
        mock_config.return_value = {'MRNA_PATH': '/mock/path'}
        
        with patch('pathlib.Path.glob') as mock_glob:
            mock_glob.return_value = []  # No files found
            
            gene = Gene.from_file('NONEXISTENT', 'hg38')
            assert gene is None