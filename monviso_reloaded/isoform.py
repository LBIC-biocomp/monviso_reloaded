
from pathlib import Path
from typing import Union
import subprocess
from Bio.Blast import NCBIWWW as blastq
from Bio.Blast import NCBIXML as blastparser
from Bio import SeqIO

from .file_handler import FileHandler
from .cobalt_wrapper import Cobalt
from .template import Template

class Isoform:
    def __init__(self,gene_name:str,sequence: list[str],isoform_index: int,out_path: Union[str,Path],mutations: list):
        self.gene_name=gene_name
        self.isoform_index=isoform_index
        self.isoform_name="isoform"+str(self.isoform_index)
        self.mutations=mutations
        self.first_line=sequence[0]
        self.sequence=sequence[1:]
        self.out_path=Path(out_path,self.isoform_name)

        #just for initialization
        self.mutation_score=0 
        self.structural_score=0 
        self.selection_score=0 
        self.templates=[]
        self.modellable= True #False if all templates get excluded
        self.aligned_sequence=""
        
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
                
    def _extract_pdb_names(self, max_pdb: int) -> list :
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
        
    def _template_alignment(self,cobalt_home: Union[str,Path]) -> None:
        """Create a unique file with the templates sequences and run a cobalt MSA.
        Args:
            cobalt_home (Union[str,Path]): Home of the Cobalt program, where executables are stored.
        """
        content=">"+self.gene_name+" "+self.isoform_name+"\n"
        content+="".join(self.sequence)+"\n"
        for template in self.templates:
            content+=">"+template.pdb_name+"_"+template.pdb_chain+"\n"
            content+=template.sequence+"\n"
            
        with FileHandler() as fh:
            templates_path=Path(self.out_path,"templates_sequences.fasta")
            fh.write_file(templates_path,content)
            aligned_path=Path(self.out_path,"templates_aligned.fasta")
            if fh.check_existence(aligned_path):
                print("Cobalt output file for templates already present in folder.")
            else:
                with Cobalt() as cobalt:
                    cobalt.run(templates_path,aligned_path,cobalt_home)
                    print(f"Cobalt alignment for {self.gene_name} {self.isoform_name} templates done.")
        
    def load_templates(self,max_pdb:int,resolution_cutoff: float,cobalt_home: Union[str,Path]):
        """Load a list of N (=max_pdb) PDB ids to use as templates.
        Create a Template object for each of the pdb files. Append
        the new object to self.templates.

        Args:
            max_pdb (int): Maximum number of PDB templates to use.
            resolution_cutoff (float): Accept CryoEM and Xray structures within
            this resolution cut-off.
            cobalt_home (Union[str,Path]): Home of the Cobalt program, where executables are stored.
        """
        templates_list=self._extract_pdb_names(max_pdb)
        for pdb in templates_list:
            template=Template(pdb, self.out_path,self.gene_name,self.isoform_name,resolution_cutoff)
            if template.usable:
                self.templates.append(template)
                
        self._template_alignment(cobalt_home)
        
       
    def calculate_mutation_score(self,mappable_mutations: list) -> None:
        """Calculate and save as local attribute the score of the mutation
           function. Calculated as: mutations that can me mapped on the isoform
           divided by the number of mutations that can be mapped on at least one
           isoform of the gene. 

        Args:
            mappable_mutations (list): The list of mutations that can be mapped
                                       on at least one isoform of the gene. 
        """
        
        score=len(self.mutations)/len(mappable_mutations)
        self.mutation_score=score

    def calculate_structural_score(self) -> None:
        """Load the sequence of the templates for the aligned file.
        Add the aligned sequence as an attribute of the template.
        Calculate how many residues of the target sequence are covered
        by at least one template.
        """
        aligned_templates_path=Path(self.out_path,"templates_aligned.fasta")
        with FileHandler() as fh:
            alignment=fh.read_file(aligned_templates_path).split(">")
            
            total_number_residues = 0 # to increase when reading alignment
            modellable_residues= 0 # to be increase when reaading alignment
            
            print([t.pdb_name for t in self.templates]) ##DEBUG all template name            
            print([t.splitlines()[0] for t in alignment[1:]])
            self.aligned_sequence="".join(alignment[1].splitlines()[1:])
            

            templates_alignment=[]
            for alignment_index,aligned_template in enumerate(alignment[2:]):
                template_sequence="".join(aligned_template.splitlines()[1:])
                templates_alignment.append(template_sequence)
                self.templates[alignment_index].add_aligned_sequence(template_sequence)
                

            for residue_index,residue in enumerate(self.aligned_sequence):
                if residue!="-":
                    total_number_residues+=1
                    if sum([template[residue_index]!="-" for template in templates_alignment])>0:
                        modellable_residues+=1
                        
            
            score=modellable_residues/total_number_residues
            self.structural_score=score
    
    def filter_templates_by_sequence_identity(self,sequence_identity_cutoff:float) -> None:
        """Take the list of templates and exclude the ones with a sequence identity
        lower or equal to the cutoff.

        Args:
            sequence_identity_cutoff (float): cut-off value
        """
        filtered_list=[]
        for template in self.templates:
            if template.sequence_identity>sequence_identity_cutoff:
                filtered_list.append(template)
            else:
                print(f"{template.sequence_identity} < {sequence_identity_cutoff}")
                print(f"Template {template.pdb_name} excluded for low sequence identity")
                
        if len(filtered_list)==0:
            print("Resolution and sequence identity cutoff let to an empty list of"+
                  f" templates for {self.gene_name} {self.isoform_name}")
            self.modellable=False
        self.templates=filtered_list
        
    def calculate_selection_score(self,w1:float=10.0,w2:float=10.0) -> None:
        """Calculate and save the selection score as attribute.
        The score is calculated with:
        score= w1*structural_score + w2*mutation_score

        Args:
            w1 (float): weight of the structural score
            w2 (float): weight of the mutation score
        """
        
        self.structural_score=w1*self.structural_score+w2*self.mutation_score