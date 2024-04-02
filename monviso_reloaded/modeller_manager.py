import os
import subprocess
from pathlib import Path

from .file_handler import FileHandler


class Modeller_manager:
    def __init__(self, isoform, mutation: list, modeller_exec: str, model_cutoff: int,
                 number_of_wt:int,number_of_mut:int):
        self.isoform = isoform
        self.mutation = mutation
        self.sequence_to_model = self.isoform.aligned_sequence[:]
        self.modeller_exec = modeller_exec
        self.model_cutoff = model_cutoff
        self.logged_scores=[]
        self.num_chains=1
        self.chain_starts=[]
        self.mutation_is_modellable=True
        
        if "".join(mutation)=="WT":
            self.num_models=number_of_wt
        else:
            self.num_models=number_of_mut
    def write(self):

        print(
            f"Modelling {self.isoform.gene_name}"
            f" {self.isoform.isoform_name} "+"".join(self.mutation)
        )
        self.write_alignment()
        self.write_script()
    
    def _check_modellability(self,sequence) -> bool:
        """Sometimes the mutation can be mapped on the isoform, but it is not
        in the range of the protein that is covered by the templates. The Gene
        object that starts the modelling process, should be updated to stop the
        modelling, and write and informed report.
        
        This function checks the mutation residue number in the sequence and returns
        a boolean representing the presence of the mutated residue. 
        
        Args:
        sequences (list): A list of [residue number, residue], with one element
        for each residue of the target sequence, after adding chain breaks.
        
        Returns:
        bool: Presence of the mutated residue in the final alignment."""
        
        if self.mutation=='WT':
            return True
        
        resnum_to_mutate=int(self.mutation[1])
        resnumbers=[r[0] for r in sequence]
        self.mutation_is_modellable = resnum_to_mutate in resnumbers
        
        if not self.mutation_is_modellable:
            self.isoform.mutations_not_in_structure.append(self.mutation)
        
        return self.mutation_is_modellable

    def write_script(self):
        alignment_name = "../modeller_input_" + "".join(self.mutation) + ".dat"
        output_name = "../modeller_output_" + "".join(self.mutation) + ".dat"
        script_path = Path(
            self.isoform.out_path,
            "run_modeller_" + "".join(self.mutation) + ".py",
        )
        template_names = [
            t.pdb_name + "_" + t.pdb_chain + "_clean"
            for t in self.isoform.templates
        ]
        content = (
            """from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb

class MyModel(AutoModel):
    def special_patches(self, aln):
        # Rename both chains and renumber the residues in each
        self.rename_segments(segment_ids="""+str(["A" for  i in range(self.num_chains)])+""",
                             renumber_residues="""+str(self.chain_starts)+""")

log.verbose()
env = environ()
env.io.atom_files_directory = '../templates/'
env.io.hetatm = True
a = MyModel(env,
    alnfile =\""""
            + alignment_name
            + """\",
    knowns = ("""
            + str(template_names)
            + """),
    sequence = \""""
            + self.isoform.gene_name
            + """\",
    assess_methods=(assess.DOPE, assess.GA341))
a.starting_model= 1
a.ending_model  = """+str(self.num_models)+"""
a.make()
ok_models = filter(lambda x: x['failure'] is None, a.outputs)
toscore = 'DOPE score'
ok_models = sorted(ok_models, key=lambda k: k[toscore])
models = [m for m in ok_models[0:10]]
myout = open(\""""
            + output_name
            + """\", "w")
for m in models:
        myout.write(str(m['name']) + " (DOPE SCORE: %.3f)" % (m[toscore]))
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib')
mdl = complete_pdb(env, m['name'])
s = selection(mdl)
s.assess_dope(output='ENERGY_PROFILE NO_REPORT', file=\""""
            + self.isoform.gene_name
            + """\", normalize_profile=True, smoothing_window=15)"""
        )

        with FileHandler() as fh:
            fh.write_file(script_path, content)

    def _mutate_reside(self, mutation, sequence:str) -> bool:
        """Take a mutation in the format
        [1 letter amino acid,residue number, 1 lett. amino acid]
        and apply it to the aligned sequence to model,
        taking into account all the "-"'s.

        Returns edited sequence
        """
 
        if mutation=="WT":
            return sequence
        

        i = int(mutation[1]) - 1
        while i < len(sequence):
            actual_residue_index = i - sequence[:i].count("-")
            if (
                actual_residue_index + 1 == int(mutation[1])
            ) and sequence[i] == mutation[0]:
                sequence = (
                    sequence[:i]
                    + mutation[2]
                    + sequence[i + 1 :]
                )
                return sequence
            i += 1

        print(f"Could not apply mutation {mutation}!")
        return sequence

    def _add_number_to_sequence(self,sequence:str)->list:
        """From a sequence, return a 2D list with [resnum,resiude]"""
        numbered_seq=[]
        i=1
        for char in sequence:
            if (char=="-"):
                numbered_seq.append(["-",char])
            else:
                numbered_seq.append([i,char])
                i+=1
        return numbered_seq
    
    def _save_chain_starts(self,sequence_with_resnum:list)->None:
        """Given a sequence saved as a list [resnum, residue],
        save all the residue numbers that represent an interruption
        in the continuity of the protein chain."""

        #lists al resnums without gaps
        filtered_numbers=[r[0] for r in sequence_with_resnum if r[0]!="-"]

   
        self.chain_starts.append(filtered_numbers[0])

        for i in range(1,len(filtered_numbers)):
            if filtered_numbers[i]!=filtered_numbers[i-1]+1:
                self.chain_starts.append(filtered_numbers[i]-1)
    
    def _add_chain_breaks(self, sequences: list) -> list:
        """For alignments in which there is no coverage for self.model_cutoff+
        residues, the section without coverage is replaced by
        chain breaks.

        Args:
            sequences (list): A list of [name, sequence] lists,
            one for each sequence to be written in the modeller
            alignment input file.

        Returns:
            list: Returns the same list, with chain breaks.
        """

        # Separate protein names and aligned sequences in two objects
        names = [object[0] for object in sequences]
        aligned_seq = [object[1] for object in sequences]

        #Add object to keep trak of residue number in target seq
        num_target_seq=self._add_number_to_sequence(aligned_seq[0])
        # Calculate worst case number of residues without coverage
        max_non_covered = max([seq.count("-") for seq in aligned_seq])

        # Search for non-covered subsequences with lenght between
        # max_non_covered and model_cutoff
        while max_non_covered >= self.model_cutoff:

            # check if "-" repeated max_non_covered times
            # is present in all the aligned structures (except target seqence)
            check_presence = [
                "-" * max_non_covered in seq for seq in aligned_seq[1:]
            ]
            if sum(check_presence) == (len(aligned_seq) - 1):

                # Find position of the sequence of "-"'s in all
                # aligned sequences for each position it will
                # check if the coverage is missing in all the
                # sequences
                positions = [
                    seq.find("-" * max_non_covered) for seq in aligned_seq[1:]
                ]
                for pos in positions:
                    check_position = [
                        seq[pos : pos + max_non_covered]
                        == "-" * max_non_covered
                        for seq in aligned_seq[1:]
                    ]
                    if sum(check_position) == (len(aligned_seq) - 1):
                        # If that is the case, replace that section in
                        # all sequences (+ target sequences) with a chain
                        # break ("/").
                        num_target_seq=num_target_seq[:pos]+num_target_seq[pos + max_non_covered :]
                        for i, seq in enumerate(aligned_seq):
                            aligned_seq[i] = (
                                seq[:pos] + "/" + seq[pos + max_non_covered :])
                            self.num_chains+=1

                        max_non_covered += (
                            1  # This repeats the check for the current value
                        )
            max_non_covered -= 1

        sequences = [[names[i], aligned_seq[i]] for i in range(len(names))]
        
        #Remove "/" at the beginning or end of the sequences:
        if sequences[0][1][0]=="/":
            sequences=[[s[0],s[1][1:]] for s in sequences]
        if sequences[0][1][-1]=="/":
            sequences=[[s[0],s[1][:-1]] for s in sequences]
        #Next line is to save resnum at chain breaks as local attributes
        self._save_chain_starts(num_target_seq)
        
        #Next line is to check if mutation can still be mapped
        _ =self._check_modellability(num_target_seq)
        
        return sequences

    def write_alignment(self):

        alignment_name = "modeller_input_" + "".join(self.mutation) + ".dat"
        output_path = Path(self.isoform.out_path, alignment_name)
        # Create an object storing all names and aligned sequences
        sequences = []
        with FileHandler() as fh:
            file_path=Path(self.isoform.out_path,"templates_aligned.fasta")
            aligned_sequence_file=fh.read_file(file_path)
            aligned_sequences=aligned_sequence_file.split(">")[1:]
            sequences=[[x.split("\n")[0].split()[-1],"".join(x.split('\n')[1:])] for x in aligned_sequences]

        #Add mutation if needed
        sequences[0][1]=self._mutate_reside(self.mutation,sequences[0][1])
        
        # Add chain breaks in place of long seqs with no coverage
        sequences = self._add_chain_breaks(sequences)


        # Start writing the content string to be printed in the file
        content = ""
        content += ">P1;" + self.isoform.gene_name + "\n"
        content += "sequence:" + self.isoform.gene_name + ":.:.:.:.::::\n"
        content += sequences[0][1] + "*\n"

        # Add templates
        for template_sequence in sequences[1:]:
            content += ">P1;" + template_sequence[0] + "_clean" + "\n"
            content += (
                "structureX:"
                + template_sequence[0]
                + "_clean"
                + ":.:.:.:.::::\n"
            )
            content += template_sequence[1] + "*\n"

        #After the last alignment it is know if the mutation
        #Can be certainly applied. Apply last check before
        #Writing
        if self.mutation_is_modellable:
            with FileHandler() as fh:
                fh.write_file(output_path, content)

    def run(self) -> None:
        
        if self.mutation_is_modellable==False:
            print("Modelling of mutation "+"".join(self.mutation)+
                  " was eexcluded on "+self.isoform.gene_name+" "+
                  self.isoform.isoform_name+" due to missing coverage.")
            self.load_log_file()
            return None
        
        model_path = Path(
            self.isoform.out_path,
            self.isoform.gene_name + "_" + "".join(self.mutation) + "_model",
        )
        with FileHandler() as fh:
            fh.create_directory(model_path)
        home_working_directory = os.getcwd()
        os.chdir(str(model_path))
        script_path = "../run_modeller_" + "".join(self.mutation) + ".py"
        command = f"{self.modeller_exec} {str(script_path)}"
        subprocess.run(
            command, shell=True, universal_newlines=True, check=True
        )
        os.chdir(home_working_directory)
        self.load_log_file()

    def load_log_file(self):
        """ Open the log file after the run of modeller.
        Look for the table with the DOPE scores.
        Save it as attribute (list).
        
        But if the mutation is not modellable, due to missing template
        converage, just save a warning string as attribute.
        """
        if self.mutation_is_modellable:
            with FileHandler() as fh:
                log_path = Path(self.isoform.out_path,"run_modeller_" + "".join(self.mutation) + ".log")
                logs=fh.read_file(log_path).splitlines()
                
                table_start_index=logs.index("Filename                          molpdf     DOPE score    GA341 score")
                table_end_index=table_start_index+logs[table_start_index:].index("")
                
                self.logged_scores=logs[table_start_index:table_end_index]
        else:
            self.logged_scores=["! The mutation was not modelled, due to missing templates for that region."]