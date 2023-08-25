import contextlib
from typing import Union
from typing import List
from typing import Literal

import os
import logging
import shutil

import requests

from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from Bio.PDB import Selection
from Bio import SeqIO

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
file_handler = logging.FileHandler("LOGFILE.log")
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

console_logger =  logging.getLogger(__name__)
console_logger.setLevel(logging.DEBUG)
console_logger.basicConfig("%(asctime)s - %(name)s - %(levelname)s - %(message)s", level=logging.DEBUG)

def warning_message(msg: str) -> None:
    """Warning message

    :param msg: Message to communicate
    """
    logger.warning(msg)
    console_logger.warning(msg)

def info_message(msg: str) -> None:
    """Information message

    :param msg: Message to communicate
    """
    logger.info(msg)
    console_logger.info(msg)


def parse_input(mutation_file: str) -> Union[List, List]:
    """Parse the list of mutations and genes from the mutation_list file.

    :param mutation_list: path to the file containing the list of mutations and genes
    :return: The list of gene and mutations, the list of genes
    """
    with open(mutation_file, "r") as my_file:
        content = my_file.read()

    blocks = [block.splitlines() for block in content.split("\n\n")]
    protein_list = []
    block = 0

    while block < len(blocks):
        protein_list.append(blocks[block][0].upper())
        block += 1

    return blocks, protein_list

def make_dir(dirname) -> None:

    with contextlib.suppress(FileExistsError):
        os.mkdir(dirname)


def make_gene_directories(blocks: List, master_directory: str) -> None:
    """Create the necessary directories
    
    :param blocks: List containing the gene names and mutations
    :param master_directory: path to the main directory
    :return:
    """
    block = 0

    while block < len(blocks):
        gene_name = blocks[block][0].upper()
        if os.path.exists(gene_name):
            message = f"Folder {gene_name} exists. Removing it and recreating."
            info_message(message)
        make_dir(gene_name)
        write_mutations_file(master_directory, gene_name, blocks, block)
        block += 1


def write_mutations_file(master_directory: str, gene_name:str, blocks: list,  block: int) -> None:
    """Write the mutations file for each gene

    :param master_directory: Main directory path
    :param gene_name: Name of the gene taken into account
    :param blocks: List of all mutations and genes
    :param block: Position on the list
    """
    folder_path = os.path.join(master_directory, gene_name)

    with open(os.path.join(folder_path, "mutations.txt"), "w") as newfile:
        mutation = 1
        while mutation < len(blocks[block]):
            newfile.write(str(blocks[block][mutation]) + "\n")
            mutation += 1


def sort_list(sub_list: List) -> List:
    """Sort a list according to the second element 
    
    i.e. a list that list the mutations has the res number as second element of the sublist
    :param sub_list:
    """
    return sorted(sub_list, key=lambda x: int(x[1]))

def import_sequence(sequencefile: str) -> any:
    """Import a sequence from a FASTA file

    :param sequencefile: path of the file containing the FASTA
    :return: fasta object
    """
    return SeqIO.read(sequencefile, "fasta")

def download_pdb(pdbid: str) -> bool:
    """Download a pdb file

    :param pdbid: ID of the PDB file
    :return: True if downloaded, False if not
    """
    pdb_url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    response = requests.get(pdb_url, timeout=5)

    if response.status_code != 200:
        warning_message = f"I could not download {pdbid.upper()}."
        logger.warning(warning_message)
        return False

    with open(f"{pdbid}.pdb", "wb", encoding="utf-8") as pdb_file:
        pdb_file.write(response.content)
    return True

def get_structure(pdbname: str) -> any:
    """Open a PDB file

    :param pdbname: Path to the PDB file (including .pdb)
    :return: Structure object
    """
    parser_pdb = PDBParser(PERMISSIVE=1)
    structure_id = "pratchett"
    return parser_pdb.get_structure(structure_id, pdbname)


def check_empty_pdb(pdbid: str) -> int:
    """Verify that the PDB file is not empty

    :param pdbid: PDBID 
    :return: Number of residues in the PDB
    """
    structure = get_structure(f"{pdbid}.pdb")
    io = PDBIO()
    io.set_structure(structure)
    residues = Selection.unfold_entities(structure, 'R')
    return len(residues)