# Run Class

The `Run` class acts as the central coordinator within the MoNvIso software tool, orchestrating the sequence of operations necessary for the analysis and modeling of protein isoforms. This class manages the execution flow, starting from reading user inputs to the final step of modeling, ensuring that each process is executed in the correct order.

## Overview

The `Run` class is designed to streamline the execution of various tasks including loading user inputs, creating Gene and Isoform objects, and managing the execution of their respective methods such as alignment, modeling, and reporting. It serves as a process manager, coordinating between user inputs, database parsing, and the execution of modeling algorithms.

## Class Methods

### `__init__(self)`

Initializes a new instance of the Run class, setting up the necessary attributes for managing the execution flow.

- **Attributes**:
  - `args`: A list to hold command line arguments.
  - `parameters`: A list to store parameters from the parameters file.
  - `mutation_list`: A list to maintain mutations parsed from the input file.
  - `input_parser`: An instance of `InputParser` for parsing user inputs.
  - `genes`: A list to hold `Gene` instances created based on the mutation list.

### `load_input(self, argv) -> None`

Loads user input from the command line and parameters file, saving them as attributes of the class.

- **Parameters**:
  - `argv`: Command line arguments.

### `load_mutation_list(self) -> None`

Parses the list of mutations and genes from the mutation_list file and saves it as an attribute.

### `create_genes(self) -> None`

Processes the mutation list to create `Gene` instances for each gene mentioned, storing them in the class's genes list.

### `create_isoforms(self) -> None`

Loads isoforms for each gene from the Uniprot database, enriching the gene instances with isoform data.

### `run_blastp(self) -> None`

Initiates a BLASTp search for every loaded isoform to find potential sequence alignments.

### `run_cobalt(self) -> None`

Executes a COBALT run for sequence alignment of the loaded isoforms.

### `run_hmmsearch(self) -> None`

Performs an HMM search for every loaded isoform, further refining the alignment process.

### `load_templates(self) -> None`

Loads structural templates for the isoforms based on specified parameters.

### `select_isoforms(self) -> None`

Selects isoforms for modeling based on weighted criteria, including structure, mutation impact, and sequence identity.

### `start_modeller(self) -> None`

Runs the Modeller tool for the selected isoforms to generate structural models.

### `write_report(self)`

Generates a comprehensive report detailing the outcomes of the modeling process for each gene.
