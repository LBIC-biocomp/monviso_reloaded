import os
import subprocess
from pathlib import Path

from .file_handler import FileHandler


class Modeller_manager:
    def __init__(self, isoform, mutation: list, modeller_exec: str):
        self.isoform = isoform
        self.mutation = mutation
        self.sequence_to_model = self.isoform.aligned_sequence[:]
        self.modeller_exec = modeller_exec

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def write(self):

        print(
            f"Modelling {self.isoform.gene_name} "
            f" {self.isoform.isoform_name} {self.mutation}"
        )
        if self.mutation != "WT":
            success_mutation = self._mutate_reside(self.mutation)
            if not success_mutation:
                raise (RuntimeError("Could not apply mutation!"))
        self.write_alignment()
        self.write_script()

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

log.verbose()
env = environ()
env.io.atom_files_directory = '../templates/'
env.io.hetatm = True
a = automodel(env,
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
a.ending_model  = 1
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

    def _mutate_reside(self, mutation) -> bool:
        """Take a mutation in the format
        [1 letter amino acid,residue number, 1 lett. amino acid]
        and apply it to the aligned sequence to model,
        taking into account all the "-"'s.

        Returns True, if succesful, else False
        """
        i = int(mutation[1]) - 1
        while i < len(self.sequence_to_model):
            actual_residue_index = i - self.sequence_to_model[:i].count("-")
            if (
                actual_residue_index + 1 == int(mutation[1])
            ) and self.sequence_to_model[i] == mutation[0]:
                self.sequence_to_model = (
                    self.sequence_to_model[:i]
                    + mutation[2]
                    + self.sequence_to_model[i + 1 :]
                )
                return True
            i += 1

        print(f"Could not apply mutation {mutation}!")
        return False

    def _add_chain_breaks(self, sequences: list) -> list:
        """For alignments inw which there is no coverage for 7+
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

        # Calculate worst case number of residues without coverage
        max_non_covered = max([seq.count("-") for seq in aligned_seq])

        # Search for non-covered subsequences with lenght between
        # max_non_covered and 7
        while max_non_covered >= 7:

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
                        for i, seq in enumerate(aligned_seq):
                            aligned_seq[i] = (
                                seq[:pos] + "/" + seq[pos + max_non_covered :]
                            )

                        max_non_covered += (
                            1  # This repeats the check for the current value
                        )
            max_non_covered -= 1

        sequences = [[names[i], aligned_seq[i]] for i in range(len(names))]
        return sequences

    def write_alignment(self):

        alignment_name = "modeller_input_" + "".join(self.mutation) + ".dat"
        output_path = Path(self.isoform.out_path, alignment_name)
        # Create an object storing all names and aligned sequences
        sequences = [[self.isoform.gene_name, self.sequence_to_model]]
        for template in self.isoform.templates:
            sequences.append(
                [
                    template.pdb_name + "_" + template.pdb_chain,
                    template.aligned_sequence,
                ]
            )

        # Add chain breaks in place of long seqs with no coverage
        sequences = self._add_chain_breaks(sequences)

        # Start writing the content string to be printed in the file
        content = ""
        content += ">P1;" + sequences[0][0] + "\n"
        content += "sequence:" + sequences[0][0] + ":.:.:.:.::::\n"
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

        with FileHandler() as fh:
            fh.write_file(output_path, content)

    def run(self) -> None:
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
