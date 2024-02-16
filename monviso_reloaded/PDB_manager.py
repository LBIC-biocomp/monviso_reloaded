from pathlib import Path
from typing import Union

from Bio.PDB import PDBIO, PDBList, PDBParser, Select, Selection
from Bio.PDB.Polypeptide import index_to_one, three_to_index

from .file_handler import FileHandler


class ChainSelection(Select):
    def __init__(self, chain_letters, standard_atoms=True):
        self.chain_letters = chain_letters.upper()
        self.standard_atoms = standard_atoms
        self.standard_residues = [
            "ALA",
            "ARG",
            "ASN",
            "ASP",
            "CYS",
            "GLU",
            "GLN",
            "GLY",
            "HIS",
            "ILE",
            "LEU",
            "LYS",
            "MET",
            "PHE",
            "PRO",
            "SER",
            "THR",
            "TRP",
            "TYR",
            "VAL",
        ]

        self.first_model = True  # see accet_model method

    def accept_model(self, model):
        # Accept only the first model
        # The attribute is initialized as True
        # and turns to False after the first run.

        if self.first_model:
            self.first_model = False
            return True
        else:
            return False

    def accept_residue(self, residue):
        # Accept only residues that are in the list of standard amino acids
        return residue.resname.strip() in self.standard_residues

    def accept_chain(self, chain):
        # Filter the chain
        return chain.id == self.chain_letters

    def accept_atom(self, atom):
        # Filter for standard atoms if requested
        if self.standard_atoms:
            return True
        return False


class PDB_manager:
    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def downloadPDB(self, pdb: str, out_path: Union[str, Path]) -> Path:
        right_pdbname = f"{pdb}.pdb"
        wrong_pdbname = f"pdb{pdb}.ent"
        download_dir = ".pdbdownloads"

        filepath = Path(out_path, download_dir, right_pdbname)
        with FileHandler() as fh:
            if not fh.check_existence(Path(out_path, download_dir)):
                fh.create_directory(Path(out_path, download_dir))

            if fh.check_existence(filepath):
                return filepath
            else:
                pdbl = PDBList()
                filename = pdbl.retrieve_pdb_file(
                    pdb,
                    pdir=str(Path(out_path, download_dir)),
                    file_format="pdb",
                    overwrite=True,
                )

                if fh.check_existence(filename):
                    fh.move_file(
                        Path(out_path, download_dir, wrong_pdbname), filepath
                    )
                    return filepath
                else:
                    raise (FileNotFoundError(f"Could not download PDB {pdb}"))

    def extract_clean_chain(
        self,
        input_pdb_path: Union[Path, str],
        output_pdb_path: Union[Path, str],
        chain_letter: str,
        resolution_cutoff: float,
    ):
        """Take an input path, save the standard atoms
        of chain 'chain_letter' in
        a filtered new PDB file, if resolution is better
        than the parameter 'resolution'.
        The exception is with the NMR structures that will
        be included in all cases.

        Args:
            input_pdb_path (Union[Path,str]): The original PDB file path
            to be filtered
            output_pdb_path (Union[Path,str]): The path of the output PDB
            chain_letter (str): Letter of the chain to extract.

        Returns:
            resolution (float or None): The resolution of the
            X-Ray or CryoEM structure
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", str(input_pdb_path))
        if structure.header["resolution"]:
            if structure.header["resolution"] <= resolution_cutoff:
                with FileHandler() as fh:
                    io = PDBIO()
                    io.set_structure(structure)
                    if not fh.check_existence(output_pdb_path):
                        io.save(
                            str(output_pdb_path), ChainSelection(chain_letter)
                        )

                    # remove HETATM from saved file
                    saved_file = fh.read_file(output_pdb_path).splitlines()
                    saved_file = "\n".join(
                        [line for line in saved_file if "HETATM" not in line]
                    )
                    fh.write_file(output_pdb_path, saved_file)
                    return structure.header["resolution"]

        else:
            # The content of the following "if" cleans the NMR structures
            # But it's still missing an estimate of the reolution. Therefore
            # these structures will be cleaned, and included in all cases..
            if "nmr" in structure.header["structure_method"]:
                with FileHandler() as fh:
                    io = PDBIO()
                    io.set_structure(structure)
                    if not fh.check_existence(output_pdb_path):
                        io.save(
                            str(output_pdb_path), ChainSelection(chain_letter)
                        )
                print(f"NMR structure in path {output_pdb_path}")
                return 0

        print(
            f"The file {str(input_pdb_path)} was "
            + "exluded due to poor resolution."
        )
        return 9999

    def extract_fasta(
        self,
        pdb_name: str,
        pdb_path: Union[str, Path],
        output_fasta_path: Union[str, Path],
    ) -> None:
        """Load a PDB file and return the fasta sequence as a string.

        Args:
            pdb_name: Name of the PDB file and chain letter
            pdb_path (Union[str,Path]): Path to the PDB structure.
            output_fasta_path (Union[str,Path]): Path where to save the
            fasta sequence.
        """
        with FileHandler() as fh:
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", str(pdb_path))
            io = PDBIO()
            io.set_structure(structure)

            residues = Selection.unfold_entities(structure, "R")
            resnames = [x.get_resname() for x in residues]
            sequence = "".join(
                [index_to_one(three_to_index(resname)) for resname in resnames]
            )
            content = ">" + pdb_name + "\n"
            content += sequence
            content += "\n"
            # This next check could be moved at the beginning
            # and if .fasta file already exists, it should
            # be loaded
            if not fh.check_existence:
                fh.write_file(output_fasta_path, content)
            return sequence
