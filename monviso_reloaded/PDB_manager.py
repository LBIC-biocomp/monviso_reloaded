from pathlib import Path
from typing import Union

from Bio.PDB import PDBIO, PDBList, PDBParser, Select, Selection, MMCIFParser
from Bio.PDB.Polypeptide import index_to_one, three_to_index

from .file_handler import FileHandler


class ChainSelection(Select):
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters.upper()
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
        
        self.standard_atoms=['C',
                             'CA',
                             'CB',
                             'CD',
                             'CD1',
                             'CD2',
                             'CE',
                             'CE1',
                             'CE2',
                             'CG',
                             'CG1',
                             'CG2',
                             'CZ',
                             'N',
                             'ND1',
                             'ND2',
                             'NE',
                             'NE2',
                             'NH1',
                             'NH2',
                             'NZ',
                             'O',
                             'OD1',
                             'OD2',
                             'OE1',
                             'OE2',
                             'OG',
                             'OG1',
                             'OH',
                             'SD',
                             'SG']

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
        return (chain.id.upper() == self.chain_letters)

    def accept_atom(self, atom):
        # Filter for standard atoms
        return atom.name.strip() in self.standard_atoms


class PDB_manager:
    def __init__(self):
        pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def downloadPDB(self, pdb: str, out_path: Union[str, Path]) -> Path:
        right_pdbname = f"{pdb}.cif"
        download_dir = ".pdbdownloads"

        filepath = Path(out_path, download_dir, right_pdbname)
        with FileHandler() as fh:
            if not fh.check_existence(Path(out_path, download_dir)):
                fh.create_directory(Path(out_path, download_dir))

            if fh.check_existence(filepath):
                return filepath
            else:
                pdbl = PDBList(server='https://files.wwpdb.org')
                filename = pdbl.retrieve_pdb_file(
                    pdb,
                    pdir=str(Path(out_path, download_dir)),
                    overwrite=True,
                    obsolete=False
                )

                if fh.check_existence(filename):
                    fh.move_file(
                        filename, filepath
                    )
                    return filepath
                else:
                    RuntimeWarning(f"Could not download structure with ID {pdb}.")

                return False

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
        if str(input_pdb_path).endswith("pdb"):
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure("structure", str(input_pdb_path))
            resolution = structure.header.get("resolution", None)
        else:  # Assume CIF format
            parser = MMCIFParser(QUIET=True)
            structure = parser.get_structure("structure", str(input_pdb_path))
            resolution = None
            if "_refine.ls_d_res_high" in parser._mmcif_dict:
                try:
                    resolution = float(parser._mmcif_dict["_refine.ls_d_res_high"][0])
                except ValueError:
                    pass  # Resolution could not be converted to a float
            if "_em_3d_reconstruction.resolution" in parser._mmcif_dict:
                try:
                    resolution = float(parser._mmcif_dict["_em_3d_reconstruction.resolution"][0])
                except ValueError:
                    pass  # Resolution could not be converted to a float
        if resolution:
            if resolution <= resolution_cutoff:
                with FileHandler() as fh:
                    io = PDBIO()
                    io.set_structure(structure)
                    if not fh.check_existence(output_pdb_path):
                        io.save(
                            str(output_pdb_path), ChainSelection(chain_letter)
                        )

                    # remove HETATM from saved file
                    #saved_file = fh.read_file(output_pdb_path).splitlines()
                    #saved_file = "\n".join(
                    #    [line for line in saved_file if "HETATM" not in line]
                    #)
                    #fh.write_file(output_pdb_path, saved_file)
            return resolution

        else:
            # The content of the following "if" cleans the NMR structures
            # But it's still missing an estimate of the reolution. Therefore
            # these structures will be cleaned, and included in all cases..
            if "NMR" in structure.header["structure_method"].upper():
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
            + "excluded due to poor resolution."
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

    def change_chain_and_save(self, input_path, chain_name, output_path):
            """
            Load a PDB file, change all chain identifiers to the specified chain_name, and save it.

            Parameters:
            - input_path: Path to the input PDB file.
            - chain_name: Single-letter string for the new chain identifier.
            - output_path: Path to save the modified PDB file.
            """
            # Load the structure
            parser = PDBParser(QUIET=True)
            io = PDBIO()
            structure = parser.get_structure('structure', input_path)

            # Modify the chain identifiers
            for model in structure:
                for chain in model:
                    chain.id = chain_name  # Change chain identifier

            io.set_structure(structure)
            io.save(output_path, ChainSelection(chain_name))
