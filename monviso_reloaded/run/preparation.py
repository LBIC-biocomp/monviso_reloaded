import argparse
import contextlib
from pathlib import Path


def parse_isoform(infile) -> list[list[str]]:
    with open(infile, "r") as f:
        content = f.read()

    return [
        block.splitlines() for block in content.split(">") if "sp|" in block
    ]


def write_canonical_isoform_file(
    path: Path, isoforms: list, counter: int
) -> None:
    """Write the non canonical FASTA sequences

    :param path: Path to the GENE folder
    :param isoforms: List of isoforms
    :param counter: Counter to keep track of the list position
    """
    with Path(path, "isoform0.fasta").open("w", encoding="utf-8") as filename:
        filename.write(f">{str(isoforms[counter][0])}" + "\n")
        for line in isoforms[counter][1:]:
            filename.write(line)
            if line != isoforms[counter][-1]:
                filename.write("\n")


def write_non_canonical_isoform_file(
    path: Path, isoforms: list, counter: int, words: str
) -> None:
    """Write the non canonical FASTA sequences

    :param path: Path to the GENE folder
    :param isoforms: List of isoforms
    :param counter: Counter to keep track of the list position
    :param words: ID of the isoform
    """
    with Path(path, f"isoform{str(words[2])}.fasta").open(
        "w", encoding="utf-8"
    ) as filename:
        filename.write(f">{str(isoforms[counter][0])}" + "\n")
        for line in isoforms[counter][1:]:
            filename.write(line)
            if line != isoforms[counter][-1]:
                filename.write("\n")


def get_isoforms_from_db(
    protein_list: list, outputs: argparse.Namespace, db_location: str
) -> None:
    """Extract the FASTA sequences from the DBs

    :param protein_list: List of gene names
    :param outputs: Path to output folder
    :param db_location: Path to DBs location
    """
    # the first part deals with the canonical isoform which is located in uniprot_sprot.fasta
    canonical_db = Path(db_location, "uniprot_sprot.fasta")
    isoforms_db = Path(db_location, "uniprot_sprot_varsplic.fasta")
    isoforms = parse_isoform(canonical_db)

    for gene in protein_list:
        name = f"GN={str(gene.upper())} "
        counter = 0
        while counter < len(isoforms):
            if name in str(isoforms[counter][0]) and "Homo" in str(
                isoforms[counter][0]
            ):
                path = Path(outputs, str(gene))
                filename = f"{str(path)}isoform0.fasta"
                write_canonical_isoform_file(path, isoforms, counter)
            counter += 1

    # the second part deals with all the other isoforms which are in the file called uniprot_sprot_varsplic.fasta
    isoforms = parse_isoform(isoforms_db)

    for gene in protein_list:
        name = f"GN={str(gene.upper())}"
        numb = len(name)
        counter = 0
        while counter < len(isoforms):
            if (
                isoforms[counter][0][-numb:] == name
                and "Homo" in isoforms[counter][0]
            ):
                path = Path(outputs, str(gene))
                words = isoforms[counter][0].split(" ")
                write_non_canonical_isoform_file(
                    path, isoforms, counter, words
                )

            counter += 1


def build_master_isoform_file(path: Path) -> None:
    """Build the file to keep track of the isoforms

    :param path: Path to the GENE directory
    """

    isolist = [file for file in path.iterdir() if "fasta" in file.name]
    isolist.sort()
    numbers = []
    for i in isolist:
        if i.name[:7] == "isoform":
            try:
                numbers.append(int(i.name[7:-6]))
            except Exception:
                numbers.append(i.name[7:-6])

    with contextlib.suppress(Exception):
        numbers.sort()
    with Path("master_isoform.txt").open("w", encoding="utf-8") as miso:
        for i in numbers:
            miso.write(f"isoform{str(i)}.fasta\n")
