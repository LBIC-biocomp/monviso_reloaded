# flake8: noqa
import argparse
import logging
from pathlib import Path
from typing import Any, Dict, List, Union

import requests
from Bio import SeqIO
from Bio.PDB import PDBIO, PDBParser, Selection

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler = logging.FileHandler("LOGFILE.log")
file_handler.setLevel(logging.DEBUG)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)

console_logger = logging.StreamHandler()
console_logger.setFormatter(formatter)
logger.addHandler(console_logger)


def log_warning_message(msg: str) -> None:
    """Log warning message

    :param msg: Message to communicate
    """
    logger.warning(msg)


def log_info_message(msg: str) -> None:
    """Log information message

    :param msg: Message to communicate
    """
    logger.info(msg)


def parse_input(mutation_file_path: argparse.Namespace) -> Union[List, List]:
    """Parse the list of mutations and genes from the mutation_list file.

    :param mutation_list: path to the file containing the list of mutations and genes
    :return: The list of gene and mutations, the list of genes
    """
    with Path(mutation_file_path).open() as my_file:
        content = my_file.read()

    blocks = [block.splitlines() for block in content.split("\n\n")]
    protein_list = []
    block = 0

    while block < len(blocks):
        protein_list.append(blocks[block][0].upper())
        block += 1

    return blocks, protein_list


def make_dir(dirname) -> None:
    """Create a directory

    :param dirname: Directory path
    """
    Path(dirname).mkdir(exist_ok=True)


def get_parameters(
    parameters_path: argparse.Namespace = "parameters.dat",
) -> Dict:
    """Collect the parameters from the parameters file if provided.

    :param parameters_path: Path to the parameters file, defaults to "parameters.dat"
    :return: Dict of keywords with the associated parameters
    """
    keywords = {
        "RESOLUTION": None,
        "SEQID": None,
        "PDB_TO_USE": None,
        "DB_LOCATION": None,
        "HMMER_HOME": None,
        "COBALT_HOME": None,
        "MODEL_CUTOFF": None,
        "NUM_OF_MOD_WT": None,
        "NUM_OF_MOD_MUT": None,
    }
    keys = list(keywords)
    if not Path(parameters_path).exists():

        error_message = f"Parameters file not found in {parameters_path}, please check the path provided."
        raise TypeError(error_message)

    with Path(parameters_path).open() as my_file:
        lines = my_file.readlines()

    for key in keys:
        for line in lines:
            if key in line:
                value = line[line.find("=") + 1 :].strip()
                keywords[key] = value
    return keywords


def make_gene_directories(
    blocks: List, output_directory: argparse.Namespace
) -> None:
    """Create the necessary directories

    :param blocks: List containing the gene names and mutations
    :param master_directory: path to the output directory
    :return:
    """
    block = 0

    while block < len(blocks):
        gene_name = blocks[block][0].upper()
        folder_path = Path(output_directory, gene_name)
        if folder_path.is_dir():
            message = f"Folder {folder_path} exists."
            log_info_message(message)
        make_dir(folder_path)
        write_mutations_file(folder_path, blocks, block)
        block += 1


def write_mutations_file(folder_path: str, blocks: list, block: int) -> None:
    """Write the mutations file for each gene

    :param output_directory: Main directory path
    :param gene_name: Name of the gene taken into account
    :param blocks: List of all mutations and genes
    :param block: Position on the list
    """
    file_path = folder_path / "mutations.txt"
    with file_path.open("w", encoding="utf-8") as newfile:
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


def import_sequence(sequencefile: str) -> Any:
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


def get_structure(pdbname: str) -> Any:
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
    residues = Selection.unfold_entities(structure, "R")
    return len(residues)


def add_arguments(parser: argparse.ArgumentParser) -> None:

    parser.add_argument(
        "-i",
        "--input_file",
        help="Path to the file genes and mutations",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-o",
        "--out_path",
        help="Path to the output folder",
        type=str,
        required=True,
    )

    parser.add_argument(
        "-pf",
        "--par_file",
        help="Input the path to the parameters file",
        type=str,
        required=False,
    )
    parser.add_argument(
        "-db",
        "--db_home",
        help="path to the DBs containing the canonical and isoform sequences",
        type=str,
        required=False,
    )

    parser.add_argument(
        "-cobalt",
        "--cobalt_home",
        help="path to COBALT bin folder",
        type=str,
        required=False,
    )

    parser.add_argument(
        "-hmmer",
        "--hmmer_home",
        help="path to HMMER bin folder",
        type=str,
        required=False,
    )

    parameters_group = parser.add_argument_group(
        "input", "manually provides the inputs"
    )
    parameters_group.add_argument(
        "-res",
        "--resolution",
        help="worst accepted resolution",
        type=float,
        default=4.50,
        required=False,
    )
    parameters_group.add_argument(
        "-seqid",
        "--sequence_identity",
        help="minimum sequence identity accepted",
        type=float,
        default=25,
        required=False,
    )
    parameters_group.add_argument(
        "-max_pdb",
        "--max_pdb_templates",
        help="maximum number of PDBs to use as templates",
        type=int,
        default=10,
        required=False,
    )

    parameters_group.add_argument(
        "-wt_mod",
        "--max_model_wt",
        help="maximum number of wild-type models to build",
        type=int,
        default=1,
        required=False,
    )

    parameters_group.add_argument(
        "-mut_mod",
        "--max_model_mut",
        help="maximum number of mutant models to build",
        type=int,
        default=1,
        required=False,
    )

    parameters_group.add_argument(
        "-rc",
        "--res_cutoff",
        help="maximum number of not crystalised residues accepted before cutting a model",
        type=int,
        default=5,
        required=False,
    )


def check_arguments(args: argparse.Namespace) -> None:
    """Verify that the necessary parameters have been provided

    :param args: passed arguments
    """
    if not args.par_file and (
        not args.db_home or not args.cobalt_home or not args.hmmer_home
    ):
        argument_error(
            "Specify a parameters file or insert the paths to the DBs, COBALT and HMMER manually"
        )
    elif args.par_file and (
        args.db_home or args.cobalt_home or args.hmmer_home
    ):
        argument_error(
            "Either specify a parameters file or insert parameters manually"
        )


def argument_error(msg: str) -> None:
    """Log error and raise it

    :param msg: Message to report
    """
    log_info_message(msg)
    raise TypeError(msg)


def print_parameters(args: argparse.Namespace, parameters: Dict) -> None:
    """Print the parameters provided

    :param parameters: Provided parameters
    """
    param = (
        f"\nResolution: {parameters['RESOLUTION']}\n"
        f"SEQ ID: {parameters['SEQID']}\n"
        f"WT MODELS TO PREPARE: {parameters['NUM_OF_MOD_WT']}\n"
        f"MUTANTS MODEL TO PREPARE: {parameters['NUM_OF_MOD_MUT']}\n"
        f"MAX PDBS AS TEMPLATES: {parameters['PDB_TO_USE']}\n"
        f"RESIDUES CUTOFF: {parameters['MODEL_CUTOFF']}\n"
        f"INPUT FILE: {args.input_file}\n"
        f"OUTPUT DIRECTORY: {args.out_path}\n"
        f"DATABASES DIRECTORY: {parameters['DB_LOCATION']}\n"
        f"COBALT: {parameters['COBALT_HOME']}\n"
        f"HMMER: {parameters['HMMER_HOME']}"
    )
    log_info_message(param)


def merge_parameters(args: argparse.Namespace) -> Dict:
    """Make the parameters dictionary if they are passed as arguments

    :param args: passed arguments
    """
    return {
        "RESOLUTION": args.resolution,
        "SEQID": args.sequence_identity,
        "PDB_TO_USE": args.max_pdb_templates,
        "DB_LOCATION": args.db_home,
        "HMMER_HOME": args.hmmer_home,
        "COBALT_HOME": args.cobalt_home,
        "MODEL_CUTOFF": args.res_cutoff,
        "NUM_OF_MOD_WT": args.max_model_wt,
        "NUM_OF_MOD_MUT": args.max_model_mut,
    }
