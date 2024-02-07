from typing import Union
from pathlib import Path
from .file_handler import FileHandler
from .PDB_manager import PDB_manager

class Template:
    def __init__(self, pdb_name: str, out_path: Union[str,Path], gene_name: str,isoform_name:str,resolution_cutoff):
        """Assign PDB attributes to this object, such as name, chain, gene_name, etc.
        Within the out_path directory of the isoform, create the templates directory,
        download the pdb file, clean it, and extract resolution and fasta sequence.

        Args:
            pdb_name (str): The pdb name followed by the chain letter (e.g., "6J8G_A").
            out_path (Union[str,Path]): Path to the isoform directory. In the way the
                            code is written, you also need access to the parent of the
                            parent.
            gene_name (str): name of the gene.
            isoform_name (str): name of the isoform.
        """
        
        self.usable=True #or False after checks
        self.pdb_name=pdb_name[:4]
        self.pdb_chain=pdb_name[-1].upper()
        self.gene_name=gene_name
        self.out_path= Path(out_path)
        self.isoform_name=isoform_name
        self.templates_directory=Path(self.out_path,"templates")
        self.pdb_filename=Path(self.templates_directory,self.pdb_name+".pdb")
        self.clean_pdb_filename=Path(self.templates_directory,pdb_name+"_clean.pdb")
        self.clean_fasta_file=Path(self.templates_directory,pdb_name+".fasta")
        self.resolution_cutoff=resolution_cutoff
        self.resolution=9999  #Overwritten with Xray and CryoEM resolution
        self.sequence=''
        
        self.get_pdb_file()
        self.get_clean_pdb_chain()
        if self.usable:
            self.get_fasta()
        
    def get_pdb_file(self) -> None:
        """Create the template directory for the PDB files if does not exist.
        Download the PDB files that could be used as templates.

        Args:
            templates_list (str): list of pdbs that could be used for modelling
        """
        with FileHandler() as fh:
            if not fh.check_existence(self.templates_directory):
                fh.create_directory(self.templates_directory)
                
            with PDB_manager() as pm:
                if not fh.check_existence(self.pdb_filename):
                    file=pm.downloadPDB(self.pdb_name,self.out_path.parent.parent)
                    fh.copy_file(file,self.templates_directory)
                    
    def get_clean_pdb_chain(self) -> None:
        with PDB_manager() as pm:
            self.resolution=pm.extract_clean_chain(self.pdb_filename,self.clean_pdb_filename,self.pdb_chain, self.resolution_cutoff)
            if self.resolution>self.resolution_cutoff:
                self.usable=False
    def get_fasta(self) -> None:
        with PDB_manager() as pm:
            self.sequence=pm.extract_fasta(self.pdb_name, self.clean_pdb_filename, self.clean_fasta_file)
            