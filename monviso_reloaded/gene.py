from pathlib import Path
from typing import Union
import subprocess
from Bio.Blast import NCBIWWW as blastq
from Bio.Blast import NCBIXML as blastparser
from Bio import SeqIO
from time import sleep

from .database_parser import DatabaseParser
from .file_handler import FileHandler
from .cobalt_wrapper import Cobalt
from .PDB_manager import PDB_manager


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
        self.out_path=Path(out_path,self.isoform_name)
        self.create_directory()
        self.save_fasta_sequence()
        
    def create_directory(self) -> None:
        """Create an empty subdirectory with the isoform index
        at the output path specified by the user, and within
        the gene folder."""

        with FileHandler() as fh:
            fh.create_directory(Path(self.out_path))

    def save_fasta_sequence(self) -> None:
        
        text_output="\n".join([">"+self.first_line]+self.sequence)
        with FileHandler() as fh:
            fh.write_file(Path(self.out_path,self.isoform_name+".fasta"),text_output)
            
    def blastp_search(self) -> None:
        """Use the isoform.fasta file saved in the directory
           to start a Blastp search. If the file already exists,
           nothing is done.
        """
        print(f"Looking for homologues of {self.gene_name} {self.isoform_name}")
        file_path=Path(self.out_path,self.isoform_name+".fasta")
        out_path=Path(self.out_path,self.isoform_name+"_hits.fasta")
        
        with FileHandler() as fh:
            if fh.check_existence(out_path):
                print(f"Blastp search output file already present in folder.")
                
            else:
                fasta_file= SeqIO.read(
                            file_path, "fasta"
                            )

                results = blastq.qblast("blastp", "swissprot", fasta_file.seq, alignments=500, word_size=6)
                blastRecord = blastparser.read(results)
                text_output=">"+fasta_file.id+"\n"+str(fasta_file.seq).strip() +"\n"
                for alignment in blastRecord.alignments:
                    for hsp in alignment.hsps:
                        text_output+=f">{alignment.hit_id}\n"
                        text_output+=str(hsp.sbjct).replace("-", "") + "\n\n"
                fh.write_file(str(file_path).replace(".fasta","_hits.fasta"),text_output)
        print("Done")
    
    def create_MSA(self,cobalt_home: Union[str,Path]) -> None :
        """Take the blastp results saved in the _hits.fasta file and use them
        as query for cobalt, if the MSA are not already present in the isoform
        folder.

        Args:
            cobalt_home (Union[str,Path]): Home of the Cobalt program, where executables are stored.
        """
        hits_path=Path(self.out_path,self.isoform_name+"_hits.fasta")
        aligned_path=Path(self.out_path,"aligned.fasta")
        with FileHandler() as fh:
            if fh.check_existence(aligned_path):
                print("Cobalt output file already present in folder.")
            else:
                with Cobalt() as cobalt:
                    cobalt.run(hits_path,aligned_path,cobalt_home)
            
    def buildHMM(self,hmmer_home: Union[str,Path]) -> Union[str,Path]:
        """Take the aligned cobalt output from the aligned.fast file
        and use it as query for  hmmbuild.

        Args:
            hmmer_home (Union[str,Path]): Home of the HMMer program, where executables are stored.
            
        Returns:
            output_path (Union[str,Path]): Path to .hmm file
        """
        print(f"Building HMM for gene {self.gene_name} {self.isoform_name}")
        output_path=Path(self.out_path,self.isoform_name+".hmm")
        aligned_path=Path(self.out_path,"aligned.fasta")
        with FileHandler() as fh:
            if fh.check_existence(output_path):
                print("HMMsearch output file already present in folder.")
                return output_path
            else:
                command = f"{hmmer_home}hmmbuild {output_path} {aligned_path}"
                subprocess.run(command, shell=True, universal_newlines=True, check=True)
                return output_path
        print("Done")

    def HMMsearch(self,hmmer_home: Union[str,Path]) -> None:
        """Take the aligned cobalt output from the aligned.fast file, build .hmm file,
        and use it as query for a hmmsearch.

        Args:
            hmmer_home (Union[str,Path]): Home of the HMMer program, where executables are stored.
        """
        with FileHandler() as fh:
            hmm_path=self.buildHMM(hmmer_home=hmmer_home)
            if not fh.check_existence(hmm_path):
                raise(FileNotFoundError(".hmm file not found."))
            
            print(f"Looking for templates for {self.gene_name} {self.isoform_name}")
            templates_path=Path(self.out_path,"possible_templates.xml")
            
            if not fh.check_existence(templates_path):
                command = f"curl -L -H 'Expect:' -H 'Accept:text/xml' -F seqdb=pdb -F\
                    seq='<{str(hmm_path)}' https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch"
                try:
                    
                    result = subprocess.run(command, shell=True, universal_newlines=True, capture_output=True, check=True)
                    output = result.stdout  # Captured output as a string
                    fh.write_file(templates_path,output)
                except subprocess.CalledProcessError as e:
                    print(f"Error executing curl command: {e}")
            else:
                print(f"Templates file for {self.gene_name} {self.isoform_name} esists.\
                    Skipping hmmsearch.")
                
    def extract_pdb_names(self, max_pdb: int) -> list :
        """Extact PDB names from the hmmsearch file found in the isoform directory,
        with the name possible templates.xml.

        Args:
            max_pdb: The maximum number of PDB templates to use to model the isoform.

        Returns:
            pdb_list (list(str)): A list of the PDBids and chains that can be used
            for the modeling.
        """
        with FileHandler() as fh:
            pdb_list_path= Path(self.out_path,"all_pdbs.dat")
            top_templates_path=Path(self.out_path,"top_templates.dat")
            templates_path=Path(self.out_path,"possible_templates.xml")
            if not (fh.check_existence(pdb_list_path) and fh.check_existence(top_templates_path)):
                templates_content=fh.read_file(templates_path).splitlines()
                pdb_list=[line[16:22] for line in templates_content if "hits name" in line ]
                templates_list=pdb_list[:max_pdb]
                fh.write_file(top_templates_path,"\n".join(templates_list))
                fh.write_file(pdb_list_path,"\n".join(pdb_list))
            templates_list=fh.read_file(top_templates_path).splitlines()
            return templates_list
        
    def get_pdb_file(self,templates_list: str) -> None:
        """Create the template directory for the PDB files if does not exist.
        Download the PDB files that could be used as templates.

        Args:
            templates_list (str): list of pdbs that could be used for modelling
        """
        templates_directory=Path(self.out_path,"templates")
        with FileHandler() as fh:
            if not fh.check_existence(templates_directory):
                fh.create_directory(templates_directory)
                
            with PDB_manager() as pm:
                for pdb in templates_list:
                    if not fh.check_existence(Path(templates_directory,pdb[:4]+".pdb")):
                        file=pm.downloadPDB(pdb[:4],self.out_path.parent.parent)
                        fh.copy_file(file,templates_directory)

    
    def get_clean_pdb_chain(self,templates_list: list,resolution: float) -> None:
        templates_directory=Path(self.out_path,"templates")
        with PDB_manager() as pm:
            for pdb in templates_list:
                pdb_name=pdb[:4]
                chain=pdb[-1]
                original_pdb_path=Path(templates_directory,pdb_name+".pdb")
                filtered_pdb_path=Path(templates_directory,pdb+"_clean.pdb")
                pm.extract_clean_chain(original_pdb_path,filtered_pdb_path,chain,resolution)

            
                
    def select_pdb(self,max_pdb:int,resolution: float) -> None :
        """Take care of the template selection.
        Execute the following steps:
        1. Extract pdb names from the hmmsearch output
        2. Download the pdb files
        3. ...

        Args:
            max_pdb (int): _description_
        """
        templates_list=self.extract_pdb_names(max_pdb)
        self.get_pdb_file(templates_list)
        self.get_clean_pdb_chain(templates_list,resolution)
       