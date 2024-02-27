# Isoform Class Documentation

## Overview

The `Isoform` class is a central component of the Monviso package, encapsulating comprehensive information about gene isoforms, including their sequences, mutations, and alignments. It is designed to facilitate sophisticated calculations and analyses on these sequences, making it an invaluable tool for computational biology researchers and developers. This class supports mutation assessment, sequence alignment, and the generation of structural models, among other functionalities.

## Attributes

Upon initialization, the `Isoform` class sets up several attributes that are crucial for its operation:

- `gene_name`: The name of the gene associated with this isoform.
- `isoform_index`: An integer representing the isoform's index or identifier.
- `isoform_name`: A string derived from the isoform index, prefixed with "isoform" for easy identification.
- `mutations`: A list of mutations associated with this isoform.
- `first_line`: The first line of the sequence, often containing a descriptor or header.
- `sequence`: The actual sequence of the isoform, typically represented as a list of strings.
- `out_path`: The output directory path where the isoform's data (sequences, templates, models) will be stored.
- `mutation_score`: Initialized to 0, intended to hold a score based on mutation analysis.
- `structural_score`: Initialized to 0, intended to represent a score derived from structural analysis.
- `selection_score`: Initialized to 0, a composite score that may incorporate various scoring metrics.
- `templates`: An initially empty list, intended to hold template data used in modeling.
- `modellable`: A boolean flag indicating whether the isoform can be modeled based on available templates.
- `aligned_sequence`: A string to store the sequence aligned against a reference or template.
- `clean_aligned_sequence`: A version of `aligned_sequence` with gaps in place of long non-covered regions.
- `modeller_run`: Initially `None`, intended to hold workflow manager related to running the Modeller software for structural modeling.

## Methods

- `create_directory(self)`: Creates a directory specific to the isoform within the gene's directory structure.
- `save_fasta_sequence(self)`: Saves the isoform's sequence in FASTA format to the specified output path.
- `blastp_search(self)`: Initiates a BLASTP search with the isoform sequence as the query, storing results locally.
- `create_MSA(self, cobalt_home)`: Generates a multiple sequence alignment (MSA) from BLASTP search results.
- `buildHMM(self, hmmer_home)`: Constructs an HMM file from the MSA for template searching.
- `HMMsearch(self, hmmer_home)`: Conducts a search for structural templates using the generated HMM file.
- `_extract_pdb_names(self, max_pdb)`: Extracts PDB names from the HMM search results, limited by `max_pdb`.
- `_template_alignment(self, cobalt_home)`: Creates an alignment of template sequences, incorporating the isoform sequence.
- `load_templates()`: Initializes template objects based on PDB names and chain identifiers.
- `_add_chain_breaks(self, sequences, model_cutoff)`: Inserts chain breaks in alignments to facilitate sequence identity calculation.
- `calculate_mutation_score(self, mappable_mutations)`: Calculates a score based on the ratio of mappable mutations.
- `calculate_structural_score(self, model_cutoff)`: Assesses coverage of the isoform sequence by template structures.
- `filter_templates_by_sequence_identity()`: Filters out templates with sequence identity below a specified threshold.
- `calculate_selection_score()`: Computes a selection score based on structural and mutation scores.
- `run_modeller(self, mutation, modeller_exec, model_cutoff)`: Manages the execution of Modeller for structural modeling.
