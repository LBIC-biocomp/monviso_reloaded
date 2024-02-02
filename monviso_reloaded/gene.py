from pathlib import Path
from typing import Union
import subprocess
from Bio.Blast import NCBIWWW as blastq
from Bio.Blast import NCBIXML as blastparser
from Bio import SeqIO

from .database_parser import DatabaseParser
from .file_handler import FileHandler
from .cobalt_wrapper import Cobalt


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
        
        text_output="\n".join([">"+self.first_line]+self.sequence)
        with FileHandler() as fh:
            fh.write_file(Path(self.out_path,self.isoform_name,self.isoform_name+".fasta"),text_output)
            
    def blastp_search(self) -> None:
        """Use the isoform.fasta file saved in the directory
           to start a Blastp search. If the file already exists,
           nothing is done.
        """
        print(f"Looking for homologues of {self.gene_name} {self.isoform_name}")
        file_path=Path(self.out_path,self.isoform_name,self.isoform_name+".fasta")
        
        with FileHandler() as fh:
            if fh.check_existence(file_path):
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
        hits_path=Path(self.out_path,self.isoform_name,self.isoform_name+"_hits.fasta")
        aligned_path=Path(self.out_path,self.isoform_name,"aligned.fasta")
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
        output_path=Path(self.out_path,self.isoform_name,self.isoform_name+".hmm")
        aligned_path=Path(self.out_path,self.isoform_name,"aligned.fasta")
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
            templates_path=Path(self.out_path,self.isoform_name,"possible_templates.xml")
            
            command = f"curl -L -H 'Expect:' -H 'Accept:text/xml' -F seqdb=pdb -F\
                seq='<{str(hmm_path)}' https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch"
            try:
                
                result = subprocess.run(command, shell=True, universal_newlines=True, capture_output=True, check=True)
                output = result.stdout  # Captured output as a string
                fh.write_file(templates_path,output)
            except subprocess.CalledProcessError as e:
                print(f"Error executing curl command: {e}")