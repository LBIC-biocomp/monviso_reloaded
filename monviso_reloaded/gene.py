from pathlib import Path
from copy import copy

from Bio.PDB.Polypeptide import index_to_one, three_to_index

from .database_parser import DatabaseParser
from .file_handler import FileHandler
from .isoform import Isoform


class Gene:
    def __init__(self, gene_mutation_block: list[list], out_path: str, sequenceRun=False,sequence=None):
        self.name = gene_mutation_block[0].upper()
        if not sequenceRun:
            self.mutations = self._standardize_mutations(gene_mutation_block[1:])
            self.sequence_name="isoform"
            self.sequences = []
        else:
            self.sequence_name=gene_mutation_block[1]
            self.mutations = self._standardize_mutations(gene_mutation_block[2])
            self.sequences=[sequence]
        self.mappable_mutations = []
        self.out_path = Path(out_path, self.name)
        self.create_directory()
        self.isoforms = []

        self.isoforms_to_model = []

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
        canonical_sequences = db_parser.get_canonical_isoforms(self.name)
        self.sequences = canonical_sequences
        noncanonical_sequences = db_parser.get_noncanonical_isoforms(self.name)
        self.sequences += noncanonical_sequences
                
    def load_isoforms(self, db_parser: DatabaseParser, sequenceRun=False) -> None:
        """Create an Isoform instance, for each sequence
        found in the databases, if at least a mutation can
         be mapped onto it.

        :param db_parser: instance of DatabaseParser.
        """
        if not sequenceRun:
            self.load_sequences(db_parser)
        for isoform_index, sequence in enumerate(self.sequences):

            # Make a list of all the mutations that can be directly mapped
            # onto this isoform
            modellable_mutations = []
            # This makes the multiline fasta into one-line:
            if sequenceRun:
                sequence_string=sequence
            else:
                sequence_string = "".join("".join(sequence[1:]).splitlines())
            for mutation in self.mutations:
                if self._check_presence_mutated_residue(
                    sequence_string, mutation
                ):
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
            if len(modellable_mutations) == 0:
                print(
                    "None of the mutations can be mapped on"
                    f" {self.name} isoform_{isoform_index},"
                )

 #           else:
            if sequenceRun:
                sequence=[self.name+" "+self.sequence_name,sequence]
            self.isoforms.append(
                    Isoform(
                        self.name,
                        self.sequence_name,
                        sequence,
                        isoform_index,
                        self.out_path,
                        modellable_mutations,
                    )
                )
            

    def _check_presence_mutated_residue(
        self, sequence: str, mutation: list
    ) -> bool:
        """Given a mutation, check if that specific residue is present at
        that position.

        Args:
            sequence (str): The sequence of the isoform expressed as
            a one-line string.
            mutation (list): The mutation expressed as a list, with
            one-letter residue names (e.g., ['R','899','A'])

        Returns:
            can_be_mutated (bool): True if mutation can be applied
            without changing residue number.
        """
        residue_ID=int(mutation[1]) - 1

        if len(sequence)<residue_ID+1: #if sequence too short
            return False

        can_be_mutated = sequence[residue_ID] == mutation[0]
        return can_be_mutated

    def _standardize_mutations(self, mutation_list) -> list:
        """Remove white spaces from mutations. Change residue names
        to one-letter format.
        Args:
            mutation_list (list): list of mutations as obtained
            from the input file.

        Returns:
            standard_mutation_list (list): standardized list of mutations.
        """
        standard_mutation_list = []
        standard_residues = [
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

        for i, mutation in enumerate(mutation_list):
            mutation = "".join(mutation.split())  # remove whitespaces
            first_resname = []
            number = []

            letter_index = 0
            while mutation[letter_index].isalpha():
                first_resname.append(mutation[letter_index])
                letter_index += 1

            while mutation[letter_index].isnumeric():
                number.append(mutation[letter_index])
                letter_index += 1

            second_resname = mutation[letter_index:]

            # Apply checks on resiude names
            first_resname = "".join(first_resname).upper()
            second_resname = "".join(second_resname).upper()

            if len(first_resname) == len(second_resname):
                if len(first_resname) == 3:
                    if (first_resname in standard_residues) and (
                        second_resname in standard_residues
                    ):
                        first_resname = index_to_one(
                            three_to_index(first_resname)
                        )
                        second_resname = index_to_one(
                            three_to_index(second_resname)
                        )
                    else:
                        raise (
                            ValueError(
                                f"Residue names in {mutation} do not"
                                " appear to be standard."
                            )
                        )
                else:
                    if len(first_resname) != 1:
                        raise (
                            ValueError(
                                f"Residue names in {mutation} must be in"
                                "one-letter or three-letter format."
                            )
                        )

            else:
                raise (
                    ValueError(
                        "Different residue name format"
                        f" in mutation {mutation}."
                    )
                )

            standard_mutation_list.append(
                [
                    "".join(first_resname),
                    "".join(number),
                    "".join(second_resname),
                ]
            )

        return standard_mutation_list

    def report_on_selected_isoforms(self):
        """Print on screen a report on the selected isoforms."""

        print(
            f"For the gene {self.name} the following models will be created,\n"
            "if the mutation is covered by the strctural templates:"
        )
        for modellable_isoform in self.isoforms_to_model:
            print(
                f"{modellable_isoform[0].isoform_name} "
                + "".join(modellable_isoform[1])
            )
            

    def select_isoforms(
        self, w1: float, w2: float, sequence_identity_cutoff: float,
        model_cutoff: int,cobalt_home
    ) -> None:
        """Start the calculation of the scores for all isoforms. The
        templates of each isoform are filtered, based on the sequence
        identity. The isoforms that still have templates
        are ordered by their selection score.

        Args:
            w1 (float): weight of the structural function
            w2 (float): weight of the mutation function
            sequence_identity_cutoff (float): templates with a seq. identity
                                              lower than this, are excluded.
        """
        for isoform in self.isoforms:
            isoform.calculate_score(self.mappable_mutations,model_cutoff,
                                    sequence_identity_cutoff,w1,w2,
                                    cobalt_home)

        self.isoforms = [
            isoform for isoform in self.isoforms if isoform.modellable
        ]
        if len(self.isoforms) == 0:
            print(f"No modellable isoform for gene {self.name}")
        else:
            self.isoforms.sort(key=lambda x: -x.selection_score)

            # Add wild type to list of isoforms to model
            self.isoforms_to_model.append([copy(self.isoforms[0]), "WT"])

            # Take note of mutations to model
            mutations_to_model = self.mutations[:]

            # Add mutations of best isoform
            for mutation in self.isoforms[0].mutations:
                self.isoforms_to_model.append([copy(self.isoforms[0]), mutation])
                mutations_to_model.remove(mutation)

            # Try other isoforms to model all mutations
            # if there are other isoforms to select from
            if len(self.isoforms) > 0:
                for mutation in mutations_to_model:
                    check_across_isoforms = [
                        mutation in isoform.mutations
                        for isoform in self.isoforms
                    ]
                    if sum(check_across_isoforms) > 0:
                        isoform_for_mutation = self.isoforms[
                            check_across_isoforms.index(True)
                        ]
                        self.isoforms_to_model.append(
                            [isoform_for_mutation, mutation]
                        )
                        mutations_to_model.remove(mutation)

            # Print the mutations with that cannot
            # be associated to any isoform
            if len(mutations_to_model) > 0:
                print(
                    "The following mutations will not be"
                    " modelled for gene " + self.name
                )
                print(",".join(["".join(mut) for mut in mutations_to_model]))

    def write_report(self):
        """ For every isoform and modeller run,
        append their information to a report.md file
        contained in the gene directory.
        """
        
        content="# LOG FILE\n"
        content+="- GENE NAME: "+self.name+"\n"
        content+="- PATH: "+str(self.out_path.absolute())+"\n"
        content+="- REQUESTED MUTATIONS:\n"
        
        for mutation in self.mutations:
            content+="-- "+"".join(mutation)+"\n"

        content+="- MAPPABLE MUTATIONS:\n"
        
        for mutation in self.mappable_mutations:
            content+="-- "+"".join(mutation)+"\n"
            
        content+="- ISOFORM SCORES:\n"
        
        for isoform in self.isoforms:
            content+="-- "+isoform.isoform_name+"\n"
            content+="--- Mutation:"
            content+=str(round(isoform.mutation_score,2))+" Structural:"
            content+=str(round(isoform.structural_score,2))+" Selection:"
            content+=str(round(isoform.selection_score,2))+"\n"
            content+="--- Modellable: "+str(isoform.modellable)+"\n"
            
        #Add a list of mutations excluded due missing structural converage
        non_covered_isoforms=[]
        for isoform in self.isoforms_to_model:
            if isoform[0].modeller_run.mutation_is_modellable==False:
                non_covered_isoforms.append(isoform)
        if len(non_covered_isoforms)>0:
            content+="\nThe following mutations are in a region not covered by templates:\n"
            for isoform in non_covered_isoforms:
                content+="".join(isoform[1])+"\n"
        
        content+="\n## MODELS\n"
        
        for isoform in self.isoforms_to_model:
            #isoform[0] is the object, isoform[1] is the mutation
            content+="- "+isoform[0].isoform_name+" "+"".join(isoform[1])+"\n"
            content+="- Templates: "
            
            for template in isoform[0].templates:
                content+=template.pdb_name+"_"+template.pdb_chain+" "
            
            content+="\n"
            
            content+="\n".join(isoform[0].modeller_run.logged_scores)+"\n\n"
            
        content=content.replace("B999","") #Shorten model file name
        with FileHandler() as fh:
            fh.write_file(Path(self.out_path,"report.md"),content)