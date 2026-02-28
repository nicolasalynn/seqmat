"""Tests for Gene class"""
import pickle
import sqlite3
import pytest
from pathlib import Path
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
    @patch('seqmat.lmdb_store.load_gene_from_lmdb', side_effect=ImportError)
    def test_from_file_success(self, _mock_lmdb, mock_config, tmp_path):
        """Test successful gene loading (pkl fallback when no genes_db)."""
        gene_data = {
            'gene_name': 'KRAS',
            'gene_id': 'ENSG00000133703',
            'rev': False,
            'chrm': '12',
            'transcripts': {}
        }
        (tmp_path / 'protein_coding').mkdir()
        pkl_path = tmp_path / 'protein_coding' / 'mrnas_ENSG00000133703_KRAS.pkl'
        with open(pkl_path, 'wb') as f:
            pickle.dump(gene_data, f)
        mock_config.return_value = {'MRNA_PATH': str(tmp_path)}
        gene = Gene.from_file('KRAS', 'hg38')
        assert gene is not None
        assert gene.gene_name == 'KRAS'
        assert gene.gene_id == 'ENSG00000133703'

    @patch('seqmat.sqlite_store.get_organism_config')
    @patch('seqmat.gene.get_organism_config')
    @patch('seqmat.lmdb_store.load_gene_from_lmdb', side_effect=ImportError)
    def test_from_file_success_sqlite(self, _mock_lmdb, mock_gene_config, mock_sqlite_config, tmp_path):
        """Test gene loading from genes.db when config has genes_db."""
        gene_data = {
            'gene_name': 'TP53',
            'gene_id': 'ENSG00000141510',
            'rev': False,
            'chrm': '17',
            'transcripts': {}
        }
        db_path = tmp_path / 'genes.db'
        conn = sqlite3.connect(str(db_path))
        conn.execute(
            "CREATE TABLE genes (gene_name TEXT, gene_id TEXT, biotype TEXT, data BLOB)"
        )
        conn.execute(
            "INSERT INTO genes (gene_name, gene_id, biotype, data) VALUES (?, ?, ?, ?)",
            (gene_data['gene_name'], gene_data['gene_id'], 'protein_coding', pickle.dumps(gene_data)),
        )
        conn.commit()
        conn.close()
        config = {'genes_db': str(db_path)}
        mock_gene_config.return_value = config
        mock_sqlite_config.return_value = config
        gene = Gene.from_file('TP53', 'hg38')
        assert gene is not None
        assert gene.gene_name == 'TP53'
        assert gene.gene_id == 'ENSG00000141510'
    
    @patch('seqmat.sqlite_store.get_organism_config')
    @patch('seqmat.gene.get_organism_config')
    def test_from_file_not_configured(self, mock_gene_config, mock_sqlite_config):
        """Test gene loading when organism not configured (patch both gene and sqlite_store paths)."""
        mock_gene_config.side_effect = ValueError("Organism not configured")
        mock_sqlite_config.side_effect = ValueError("Organism not configured")
        
        gene = Gene.from_file('KRAS', 'hg38')
        assert gene is None
    
    @patch('seqmat.sqlite_store.get_organism_config')
    @patch('seqmat.gene.get_organism_config')
    def test_from_file_not_found(self, mock_gene_config, mock_sqlite_config):
        """Test gene loading when file not found (no genes_db so pickle path is used, then no files)."""
        mock_gene_config.return_value = {'MRNA_PATH': '/mock/path'}
        mock_sqlite_config.return_value = {'MRNA_PATH': '/mock/path'}  # no genes_db
        
        with patch('pathlib.Path.glob') as mock_glob:
            mock_glob.return_value = []  # No files found
            
            gene = Gene.from_file('NONEXISTENT', 'hg38')
            assert gene is None