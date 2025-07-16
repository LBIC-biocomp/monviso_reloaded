import pytest
from pathlib import Path
from monviso_reloaded.input_parser import InputParser

@pytest.fixture
def temp_files(tmp_path):
    mutation_file = tmp_path / "mutations.txt"
    sequence_file = tmp_path / "sequences.txt"
    docking_file = tmp_path / "docking_partners.txt"
    parameter_file = tmp_path / "params.dat"

    mutation_file.write_text("""GRIN1
R       844     C
Ala     349     Thr
Pro     578     Arg
Ser     688     Tyr
Tyr     647     Ser

GRIN2B
E413G
C436R
M1342R
L1424F
PRO1439ALA #Comment
pro1439ala

#Comment
GRINERR
Val     123
Gly 123 Ala Extra
Ser	688	Tyr 
Tyr  647    Ser
""")

    sequence_file.write_text("""
>PROJ1:seq1:R3A,T7C
MKRISTTITTTITITTGNGAG

>PROJ2:seq2:
MPRVLLISLVTAALAAGQGDS
""")

    docking_file.write_text("""
GRIN1
GRIN2B

PROJ1
PROJ2
""")

    parameter_file.write_text("""
DB_LOCATION=/work/run
MEGADOCK_HOME=/work/Monviso/MEGADOCK
HDOCKLITE_HOME=/work/Monviso/HDOCKlite-v1.1
HADDOCK_HOME=/work/Monviso/haddock3
HADDOCK_SELECTION=sasa
HADDOCK_CUTOFF=0.9
MSMS_HOME=/work/Monviso/msms
COBALT_HOME=/work/Monviso/cobalt/bin
HMMER_HOME=/work/Monviso/hmmer/bin
PESTO_HOME=/work/Monviso/PeSTo
MODELLER_EXEC=mod10.5
RESOLUTION=4.50
SEQID=25
HMM_TO_IMPORT=100
MODEL_CUTOFF=5
PDB_TO_USE=10
NUM_OF_MOD_WT=1
NUM_OF_MOD_MUT=1
W_STRUCT=10
W_MUT=10
""")

    return mutation_file, sequence_file, docking_file, parameter_file

def test_parse_mutation_file(temp_files):
    mutation_file, _, _, _ = temp_files
    parser = InputParser()
    blocks = parser.parse_input(mutation_file, expected_block_length=None)
    print(blocks[0])
    assert len(blocks) == 3
    assert blocks[0][0] == "GRIN1"
    assert "E413G" in blocks[1]
    assert "PRO1439ALA"==blocks[1][5]

def test_parse_sequence_file(temp_files):
    _, sequence_file, _, _ = temp_files
    parser = InputParser()
    parsed = parser.parse_sequences(sequence_file)
    assert parsed[0] == ["PROJ1", "seq1", ["R3A", "T7C"]]
    assert parsed[1] == "MKRISTTITTTITITTGNGAG"
    print (parsed[2])
    assert parsed[2] == ["PROJ2", "seq2", []]
    assert parsed[3] == "MPRVLLISLVTAALAAGQGDS"

def test_parse_docking_file(temp_files):
    _, _, docking_file, _=temp_files
    parser = InputParser()
    parsed = parser.parse_input(docking_file,expected_block_length=2)

def test_get_parameters(temp_files):
    _, _, _,param_file = temp_files
    parser = InputParser()
    params = parser.get_parameters(param_file)
    assert params["HMMER_HOME"] == "/work/Monviso/hmmer/bin"
    assert params["SEQID"] == "25"
    assert params["NUM_OF_MOD_MUT"] == "1"
