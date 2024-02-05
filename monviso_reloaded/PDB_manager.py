from Bio.PDB import PDBList
from .file_handler import FileHandler
from typing import Union
from pathlib import Path

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
            