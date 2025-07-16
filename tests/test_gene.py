import pytest
from monviso_reloaded.gene import Gene

# Mock sequence and mutation data
SIMPLE_MUTATIONS = ["ALA123THR", "PHE456SER"]  # 3-letter format


BAD_MUTATIONS = ["123GLY",  # invalid
    "SER9999GLY"]

@pytest.fixture
def mock_sequence():
    # Sequence long enough to include the mutations
    return "M" * 122 + "A" + "M" * (456 - 123 - 1) + "F" + "M" * 1000

def test_standardize_mutations():
    g = Gene(["GENE1"] + SIMPLE_MUTATIONS, out_path="/tmp", sequenceRun=False)
    assert g.mutations == [["A", "123", "T"], ["F", "456", "S"]]

def test_check_presence_mutated_residue_true(mock_sequence):
    g = Gene(["GENE1"] + SIMPLE_MUTATIONS, out_path="/tmp", sequenceRun=True, sequence=mock_sequence)
    mutation = ["A", "123", "T"]
    assert g._check_presence_mutated_residue(mock_sequence, mutation)

def test_check_presence_mutated_residue_false(mock_sequence):
    g = Gene(["GENE1"] + SIMPLE_MUTATIONS, out_path="/tmp", sequenceRun=True, sequence=mock_sequence)
    mutation = ["G", "1", "A"]  # Should be M at position 1
    assert not g._check_presence_mutated_residue(mock_sequence, mutation)

def test_standardize_mutations_bad():
    with pytest.raises(ValueError):
        Gene(["GENE1"] + BAD_MUTATIONS, out_path="/tmp", sequenceRun=False)

def test_write_report_creates_file(tmp_path):
    gene = Gene(["GENE1"] + SIMPLE_MUTATIONS, out_path=tmp_path, sequenceRun=False)
    gene.mappable_mutations = [["A", "123", "T"]]
    gene.isoforms = []
    gene.isoforms_to_model = []
    gene.write_report()
    report_path = tmp_path / "GENE1" / "report.md"
    assert report_path.exists()
    content = report_path.read_text()
    assert "GENE1" in content
    assert "REQUESTED MUTATIONS" in content
