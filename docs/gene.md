# Gene Class

## Overview

The `Gene` class is a core component of the Monviso bioinformatics tool, 
designed to manage gene-specific operations. This class is responsible for 
managing the analysis of a single gene within the broader context of 
Monviso's genetic analysis capabilities. It performs tasks such as creating 
directories for gene-specific data, managing isoform objects that handle 
detailed calculations, and generating reports on the analysis.

## Attributes

- `name`: The name of the gene being analyzed.
- `mutations`: A list of mutations associated with the gene, standardized for analysis.
- `mappable_mutations`: A list of mutations that can be mapped onto at least one isoform of the gene.
- `out_path`: The output path where the gene-specific directory is created.
- `sequences`: The list of canonical and noncanonical sequences associated with the gene.
- `isoforms`: A collection of `Isoform` objects representing the different isoforms of the gene.
- `isoforms_to_model`: A list of isoforms selected for modeling based on certain criteria.

## Methods

### `__init__(self, gene_mutation_block: list[list], out_path: str)`

Initializes a new instance of the `Gene` class.

#### Parameters:

- `gene_mutation_block`: A list containing the gene name and associated mutations.
- `out_path`: The base directory where gene-specific directories will be created.

### `create_directory(self) -> None`

Creates an empty directory for the gene in the specified output path.

### `load_sequences(self, db_parser: DatabaseParser) -> None`

Loads canonical and noncanonical sequences for the gene from the Uniprot databases.

#### Parameters:

- `db_parser`: An instance of `DatabaseParser` used to fetch gene sequences.

### `load_isoforms(self, db_parser: DatabaseParser) -> None`

Creates an `Isoform` object for each sequence if at least one mutation can be mapped onto it.

#### Parameters:

- `db_parser`: An instance of `DatabaseParser`.

### `select_isoforms(self, w1: float, w2: float, sequence_identity_cutoff: float, model_cutoff: int) -> None`

Selects isoforms for modeling based on scores, sequence identity, and other criteria.

#### Parameters:

- `w1`: Weight of the structural function.
- `w2`: Weight of the mutation function.
- `sequence_identity_cutoff`: Sequence identity threshold for template filtering.
- `model_cutoff`: Cutoff value for model selection.

### `write_report(self)`

Generates and writes a report summarizing the analysis for the gene, including information about requested mutations, mappable mutations, isoform scores, and models.

## Private Methods

These methods are used internally by the `Gene` class and are not intended for direct external usage.

### `_check_presence_mutated_residue(self, sequence: str, mutation: list) -> bool`

Checks if a specific residue mutation can be applied to the sequence.

### `_standardize_mutations(self, mutation_list) -> list`

Standardizes mutation representations by removing whitespace and converting residue names to a one-letter format.

### `_report_on_selected_isoforms(self)`

Generates a console report on the isoforms selected for modeling.
