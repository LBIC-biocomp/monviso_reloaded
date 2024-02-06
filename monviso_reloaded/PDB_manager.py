from Bio.PDB import PDBList
from .file_handler import FileHandler
from typing import Union
from pathlib import Path
from Bio.PDB import PDBParser, PDBIO, Select

class ChainSelection(Select):
    def __init__(self, chain_letters, standard_atoms=True):
        self.chain_letters = chain_letters.upper()
        self.standard_atoms = standard_atoms

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
    
    def downloadPDB(self, pdb: str,out_path: Union[str,Path])-> Path:
        right_pdbname = f'{pdb}.pdb'
        wrong_pdbname = f"pdb{pdb}.ent"
        download_dir=".pdbdownloads"
        
        filepath=Path(out_path,download_dir,right_pdbname)
        with FileHandler() as fh:
            if not fh.check_existence(Path(out_path,download_dir)):
                fh.create_directory(Path(out_path,download_dir))
            
            if fh.check_existence(filepath):
                return filepath
            else:
                pdbl = PDBList()
                filename = pdbl.retrieve_pdb_file(pdb, pdir=str(Path(out_path,download_dir)),
                                                  file_format="pdb", overwrite=True)
                
                if fh.check_existence(filename):
                    fh.move_file(Path(out_path,download_dir,wrong_pdbname),filepath)
                    return filepath
                else:
                    raise(FileNotFoundError(f"Could not download PDB {pdb}"))
    
    def extract_clean_chain(self,input_pdb_path: Union[Path,str],output_pdb_path: Union[Path,str],chain_letter: str,resolution_cutoff: float):
        """Take an input path, save the standard atoms of chain 'chain_letter' in
        a filtered new PDB file, if resolution is better than the parameter 'resolution'.

        Args:
            input_pdb_path (Union[Path,str]): The original PDB file path to be filtered
            output_pdb_path (Union[Path,str]): The path of the output PDB
            chain_letter (str): Letter of the chain to extract.
        
        Returns:
            resolution (float or None): The resolution of the X-Ray or CryoEM structure
        """
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("structure", str(input_pdb_path))
        if structure.header["resolution"] <= resolution_cutoff:
            io = PDBIO()
            io.set_structure(structure)
            io.save(str(output_pdb_path), ChainSelection(chain_letter))
            return structure.header["resolution"] 
        else:
            print(f"The file {str(input_pdb_path)} was exluded due to poor resolution.")
            return None