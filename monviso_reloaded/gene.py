from typing import Any, Dict, List, Union
from .file_handler import FileHandler
from .database_parser import DatabaseParser
from pathlib import Path

class Gene:
    def __init__(self,gene_mutation_block: List[List],out_path: str):
        self.name=gene_mutation_block[0].upper()
        self.mutations=gene_mutation_block[1:]
        self.out_path= Path(out_path,self.name)
        print(self.out_path.absolute())
        self.create_directory()
        self.isoforms=[]
        
    def create_directory(self) -> None:
        """Create an empty directory with the gene name at the output path specified by the user.
        """
        with FileHandler() as fh:
            fh.create_directory(self.out_path)
    
    def load_isoforms(self,db_parser: DatabaseParser) ->None:
        print(db_parser.get_canonical_isoforms(self.name))
        print(db_parser.get_noncanonical_isoforms(self.name))

class Isoform:
    def __init__(self):
        pass