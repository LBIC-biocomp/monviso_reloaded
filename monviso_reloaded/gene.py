from pathlib import Path
from typing import Union

from .database_parser import DatabaseParser
from .file_handler import FileHandler


class Gene:
    def __init__(self, gene_mutation_block: list[list], out_path: str):
        self.name = gene_mutation_block[0].upper()
        self.mutations = gene_mutation_block[1:]
        self.out_path = Path(out_path, self.name)
        print(self.out_path.absolute())
        self.create_directory()
        self.sequences = []
        self.isoforms=[]

    def create_directory(self) -> None:
        """Create an empty directory with the gene name
        at the output path specified by the user."""

        with FileHandler() as fh:
            fh.create_directory(self.out_path)

    def load_sequences(self, db_parser: DatabaseParser) -> None:
        """Take the sequences from the Uniprot databases and
        save them as attributes

        :param db_parser: instance of DatabaseParser
        """
        canonical_sequences=db_parser.get_canonical_isoforms(self.name)
        self.sequences=canonical_sequences
        noncanonical_sequences=db_parser.get_noncanonical_isoforms(self.name)
        self.sequences+=noncanonical_sequences
        
    def load_isoforms(self,db_parser: DatabaseParser) -> None:
        """Create an Isoform instance, for each sequence found in the databases

        :param db_parser: instance of DatabaseParser.
        """
        self.load_sequences(db_parser)
        for isoform_index,sequence in enumerate(self.sequences):
            self.isoforms.append(Isoform(self.name,sequence,isoform_index,self.out_path))


class Isoform:
    def __init__(self,gene_name:str,sequence: list[str],isoform_index: int,out_path: Union[str,Path]):
        self.gene_name=gene_name
        self.isoform_index=isoform_index
        self.isoform_name="isoform"+str(self.isoform_index)
        self.first_line=sequence[0]
        self.sequence=sequence[1:]
        self.out_path=out_path
        self.create_directory()
        self.save_fasta_sequence()
        
    def create_directory(self) -> None:
        """Create an empty subdirectory with the isoform index
        at the output path specified by the user, and within
        the gene folder."""

        with FileHandler() as fh:
            fh.create_directory(Path(self.out_path,self.isoform_name))

    def save_fasta_sequence(self) -> None:
        
        text_output="\n".join([self.first_line]+self.sequence)
        with FileHandler() as fh:
            fh.write_file(Path(self.out_path,self.isoform_name,self.isoform_name+".fasta"),text_output)      
    
