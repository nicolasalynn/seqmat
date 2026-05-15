"""Tests for the (chromosome, position) -> gene lookup."""
import pickle
import sqlite3
from unittest.mock import patch

import numpy as np
import pytest

from seqmat.locator import (
    _LocationIndex,
    _coerce_pos,
    _normalize_chrm,
    gene_names_at_position,
)


def _make_gene_blob(name: str, gene_id: str, chrm: str, start: int, end: int) -> bytes:
    return pickle.dumps({
        "gene_name": name,
        "gene_id": gene_id,
        "chrm": chrm,
        "gene_start": start,
        "gene_end": end,
        "rev": False,
        "transcripts": {},
    })


@pytest.fixture
def fake_db(tmp_path):
    """A small genes.db with a handful of overlapping and adjacent genes."""
    db_path = tmp_path / "genes.db"
    conn = sqlite3.connect(str(db_path))
    conn.execute("CREATE TABLE genes (gene_name TEXT, gene_id TEXT, biotype TEXT, data BLOB)")
    rows = [
        ("KRAS", "ENSG_KRAS", "protein_coding", _make_gene_blob("KRAS", "ENSG_KRAS", "12", 100, 200)),
        ("LRMP", "ENSG_LRMP", "protein_coding", _make_gene_blob("LRMP", "ENSG_LRMP", "12", 150, 175)),
        ("FOO",  "ENSG_FOO",  "protein_coding", _make_gene_blob("FOO",  "ENSG_FOO",  "12", 300, 400)),
        ("TP53", "ENSG_TP53", "protein_coding", _make_gene_blob("TP53", "ENSG_TP53", "17", 500, 600)),
        # Unnamed gene: gene_name is empty, only gene_id available (Ensembl pseudogenes etc.)
        ("",     "ENSG_X",    "lncRNA",         _make_gene_blob("",     "ENSG_X",    "X", 1000, 2000)),
    ]
    conn.executemany(
        "INSERT INTO genes (gene_name, gene_id, biotype, data) VALUES (?, ?, ?, ?)",
        rows,
    )
    conn.commit()
    conn.close()
    return db_path


def test_normalize_chrm():
    assert _normalize_chrm("chr12") == "12"
    assert _normalize_chrm("Chr12") == "12"
    assert _normalize_chrm("CHRX") == "X"
    assert _normalize_chrm("12") == "12"
    assert _normalize_chrm(" chr1 ") == "1"


def test_coerce_pos():
    assert _coerce_pos(100) == (100, 100)
    assert _coerce_pos((100, 200)) == (100, 200)
    assert _coerce_pos([200, 100]) == (100, 200)  # reversed range
    with pytest.raises(ValueError):
        _coerce_pos((1, 2, 3))


def test_build_and_query(fake_db):
    idx = _LocationIndex.build_from_db(fake_db)
    # Point inside KRAS only
    assert idx.query_names("12", 110) == ["KRAS"]
    # Point inside both KRAS and LRMP (LRMP is nested inside KRAS in this fixture)
    assert sorted(idx.query_names("12", 160)) == ["KRAS", "LRMP"]
    # Point outside any gene
    assert idx.query_names("12", 250) == []
    # Range overlapping KRAS+LRMP and FOO
    names = idx.query_names("12", (180, 350))
    assert "KRAS" in names and "FOO" in names
    # Chrom normalization
    assert idx.query_names("chr12", 110) == ["KRAS"]
    # Unknown chrom -> empty
    assert idx.query_names("99", 100) == []
    # Empty gene_name falls back to gene_id
    assert idx.query_names("X", 1500) == ["ENSG_X"]


def test_save_and_load_roundtrip(fake_db, tmp_path):
    idx = _LocationIndex.build_from_db(fake_db)
    out = tmp_path / "gene_locations.npz"
    idx.save(out)
    assert out.exists()
    loaded = _LocationIndex.load(out)
    assert loaded.query_names("12", 110) == ["KRAS"]
    assert sorted(loaded.query_names("12", 160)) == ["KRAS", "LRMP"]


def test_gene_names_at_position_via_organism(fake_db, tmp_path):
    """End-to-end through the module-level helper, using a mocked organism config."""
    # Clear cache so a fresh build runs for this organism key
    from seqmat import locator
    locator._INDEX_CACHE.pop("hg38_test", None)
    with patch("seqmat.locator.get_organism_config") as cfg, \
         patch("seqmat.locator.get_default_organism", return_value="hg38_test"):
        cfg.return_value = {"genes_db": fake_db}
        names = gene_names_at_position("12", 110)
        assert names == ["KRAS"]
        names = gene_names_at_position("chr12", (150, 350))
        assert "KRAS" in names and "FOO" in names
