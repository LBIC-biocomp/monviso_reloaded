# Template Class

## Overview

The `Template` class is designed to manage protein templates in the context of Monviso workflow. It handles the selection, preparation, and analysis of protein database (PDB) files for modeling purposes. This includes downloading PDB files, cleaning them to focus on relevant chains, extracting sequences, via the PDB_manager class, and assessing template usability based on resolution and sequence identity.

## Methods

### `__init__(self, pdb_name: str, out_path: Union[str, Path], gene_name: str, isoform_name: str, resolution_cutoff)`
Initializes a `Template` object with specified attributes including PDB name, output path, gene name, isoform name, and resolution cutoff. It prepares the environment for template handling, including directory setup and initial file management.

### `get_pdb_file(self)`
Creates the directory for templates if it does not exist and downloads the PDB file if it's not already present. This method ensures that the required PDB files are available locally for further processing.

### `get_clean_pdb_chain(self)`
Extracts the relevant chain from the PDB file, cleaning it to remove unnecessary atoms and saving it as a new file. This process is essential for focusing on the relevant part of the protein structure.

### `get_fasta(self)`
Extracts the FASTA sequence from the cleaned PDB file of the chain of interest. This sequence is essential for further analyses, such as sequence alignment and homology modeling.

### `add_aligned_sequence(self, aligned_sequence: str)`
Saves the aligned sequence as an attribute of the `Template` object. This method is typically called after sequence alignment has been performed, linking the alignment results back to the template.

### `add_clean_aligned_sequence(self, clean_aligned_sequence: str)`
Similar to `add_aligned_sequence`, but for aligned sequences that include considerations for chain breaks or other structural nuances.

### `calculate_sequence_identity(self, reference_sequence: str)`
Calculates the sequence identity between the template's sequence and a reference sequence. This metric is crucial for evaluating the suitability of a template for modeling based on how closely it matches the target sequence.
