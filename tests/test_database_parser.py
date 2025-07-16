import pytest
import os
import datetime
from pathlib import Path
from monviso_reloaded.database_parser import DatabaseParser

FASTA_CANONICAL = """>sp|P12345|GENE1_HUMAN Some canonical description OS=Homo sapiens GN=GENE1
MSEQENCEONE
MORESEQONE

>sp|P12346|GENE2_HUMAN Another canonical OS=Homo sapiens GN=GENE2
MSEQTWO
"""

FASTA_ISOFORMS = """>sp|Q12345|GENE1VARIANT_HUMAN Variant isoform OS=Homo sapiens GN=GENE1
VSEQONE

>sp|Q12346|GENE2VARIANT_HUMAN Variant isoform OS=Homo sapiens GN=GENE2
VSEQ2
"""

@pytest.fixture
def mock_db(tmp_path):
    # Create mock canonical and isoform databases
    db_dir = tmp_path / "db"
    db_dir.mkdir()

    canonical = db_dir / "uniprot_sprot.fasta"
    isoforms = db_dir / "uniprot_sprot_varsplic.fasta"

    canonical.write_text(FASTA_CANONICAL)
    isoforms.write_text(FASTA_ISOFORMS)

    return db_dir

def test_parse_database(mock_db, capsys):
    # Use context manager for warning check
    with DatabaseParser(str(mock_db)) as db:
        assert isinstance(db.canonical_db, list)
        assert isinstance(db.isoforms_db, list)
        assert len(db.canonical_db) == 2
        assert len(db.isoforms_db) == 2

def test_get_canonical_isoforms(mock_db):
    db = DatabaseParser(str(mock_db))
    results = db.get_canonical_isoforms("GENE1")
    assert len(results) == 1
    assert "MSEQENCEONE" in results[0][1]

def test_get_noncanonical_isoforms(mock_db):
    db = DatabaseParser(str(mock_db))
    results = db.get_noncanonical_isoforms("GENE2")
    assert len(results) == 1
    assert "VSEQ2" in results[0][1]

def test_check_file_age_warns(tmp_path, capsys, monkeypatch):
    # Create fake database files
    db_path = tmp_path / "db"
    db_path.mkdir()
    fasta = db_path / "uniprot_sprot.fasta"
    varsplic = db_path / "uniprot_sprot_varsplic.fasta"
    fasta.write_text(">sp|P12345|GENE1_HUMAN OS=Homo sapiens GN=GENE1\nSEQ")
    varsplic.write_text(">sp|Q12345|GENE1VAR_HUMAN OS=Homo sapiens GN=GENE1\nSEQ")

    # Fake time 25 weeks ago
    old_time = (datetime.datetime.now() - datetime.timedelta(weeks=25)).timestamp()

    # Patch FileHandler.get_date to always return old_time
    from monviso_reloaded import file_handler
    monkeypatch.setattr(file_handler.FileHandler, "get_date", lambda self, path: old_time)

    with DatabaseParser(str(db_path)):
        pass

    captured = capsys.readouterr()
    assert "WARNING" in captured.out
    assert "older than 24 weeks" in captured.out
