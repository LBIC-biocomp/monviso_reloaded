from pathlib import Path
from Bio.PDB.Polypeptide import three_to_index,index_to_one

from .database_parser import DatabaseParser
from .file_handler import FileHandler
from .isoform import Isoform

class Gene:
    def __init__(self, gene_mutation_block: list[list], out_path: str):
        self.name = gene_mutation_block[0].upper()
        self.mutations = self._standardize_mutations(gene_mutation_block[1:])
        self.mappable_mutations=[]
        self.out_path = Path(out_path, self.name)
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
        """Create an Isoform instance, for each sequence found in the databases,
        if at least a mutation can be mapped onto it.

        :param db_parser: instance of DatabaseParser.
        """
        self.load_sequences(db_parser)
        for isoform_index,sequence in enumerate(self.sequences):
            
            # Make a list of all the mutations that can be directly mapped
            # onto this isoform
            modellable_mutations=[]
            # This makes the multiline fasta into one-line:
            sequence_string="".join("".join(sequence[1:]).splitlines()) 
            for mutation in self.mutations:
                if self._check_presence_mutated_residue(sequence_string,mutation):
                    # The mutation can be mapped, append it to the list.
                    # The list will be passed to the Isoform object.
                    modellable_mutations.append(mutation)
                    
                    # Take note of the mutations that can be mapped on
                    # at least on isoform. This will be necessary to 
                    # calculate the score of the mutation function.
                    # See doi: 10.3389/fchem.2022.1059593
                    
                    if mutation not in self.mappable_mutations:
                        self.mappable_mutations.append(mutation)
                    
            # Check if no mutations can be mapped. Skip isoform.
            if len(modellable_mutations)==0:
                print(f"None of the mutations can be mapped on {self.name} isoform_{isoform_index}")
            
            else:
                self.isoforms.append(Isoform(self.name,sequence,isoform_index,self.out_path,modellable_mutations))
            
    def _check_presence_mutated_residue(self,sequence: str, mutation: list) -> bool:
        """Given a mutation, check if that specific residue is present at
        that position.

        Args:
            sequence (str): The sequence of the isoform expressed as a one-line string.
            mutation (list): The mutation expressed as a list, with one-letter residue names
                             (e.g., ['R','899','A'])
            
        Returns:
            can_be_mutated (bool): True if mutation can be applied without changing
                residue number.
        """
        can_be_mutated= sequence[int(mutation[1])-1]==mutation[0]
        return can_be_mutated
    
    def _standardize_mutations(self,mutation_list) -> list:
        """Remove white spaces from mutations. Change residue names to one-letter format.
        Args:
            mutation_list (list): list of mutations as obtained from the input file.

        Returns:
            standard_mutation_list (list): standardized list of mutations.
        """
        standard_mutation_list=[]
        standard_residues = [
            "ALA", "ARG", "ASN", "ASP", "CYS",
            "GLU", "GLN", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO",
            "SER", "THR", "TRP", "TYR", "VAL"
        ]
        
        for i, mutation in enumerate(mutation_list):
            mutation=''.join(mutation.split()) #remove whitespaces
            
            first_resname=[]
            number=[]
            
            letter_index=0
            while mutation[letter_index].isalpha():
                first_resname.append(mutation[letter_index])
                letter_index+=1
            
            while mutation[letter_index].isnumeric():
                number.append(mutation[letter_index])
                letter_index+=1
            
            second_resname=mutation[letter_index:]
            
            # Apply checks on resiude names
            first_resname="".join(first_resname).upper()
            second_resname="".join(second_resname).upper()
            
            if (len(first_resname)==len(second_resname)):
                if len(first_resname)==3:
                    if (first_resname in standard_residues) and (second_resname in standard_residues):
                        first_resname=index_to_one(three_to_index(first_resname))
                        second_resname=index_to_one(three_to_index(second_resname))
                    else:
                        raise(ValueError(f"Residue names in {mutation} do not appear to be standard."))
                else:
                    if len(first_resname)!=1:
                         raise(ValueError(f"Residue names in {mutation} must be in one-letter or three-letter format."))  
                    
            else:
                raise(ValueError(f"Different residue name format in mutation {mutation}."))
            
            standard_mutation_list.append(["".join(first_resname),"".join(number),"".join(second_resname)])
    
        return standard_mutation_list
            
            
    def select_isoforms(self,w1: float,w2: float,sequence_identity_cutoff: float) -> None:
        for isoform in self.isoforms:

            isoform.calculate_mutation_score(self.mappable_mutations)
            isoform.calculate_structural_score()
            for template in isoform.templates:
                template.calculate_sequence_identity(isoform.aligned_sequence)
            isoform.filter_templates_by_sequence_identity(sequence_identity_cutoff)
            isoform.calculate_selection_score(w1,w2)
