# InputParser Class

The `InputParser` class extends `argparse.ArgumentParser` to provide a 
flexible and user-friendly interface for handling command line arguments as 
well as parameters from a file. It is specifically tailored for applications 
that require a dynamic configuration of execution parameters, supporting 
both direct command line input and the use of a parameters file (loaded with
the -p flag).

## Overview

This class allows for specifying computational parameters and execution 
environment settings through the command line or by providing a path to a 
structured parameters file. It is used internally by the Run class to do so.
It supports validation and merging of parameters 
to ensure that the execution environment is correctly set up according to 
the user's specifications.

## Expected Parameters

The `InputParser` class processes and validates a variety of input parameters essential for the execution of our computational workflow. Here's an overview of these parameters and their significance:

- **Resolution**: The minimum template resolution (in Å) for the modeling process, specified as `RESOLUTION`. Only valid for Cryo-EM and X-Ray structures. All NMR structures are accepted.
- **SEQ ID**: The minimum sequence identity threshold, expressed as a percentage (`SEQID`), for selecting templates. Templates with a sequence identity lower than this threshold will be discarded.
- **WT Models to Prepare**: The number of wild-type models to be prepared, indicated by `NUM_OF_MOD_WT`.
- **Mutants Model to Prepare**: Specifies the number of mutant models to prepare, `NUM_OF_MOD_MUT`.
- **Max PDBs as Templates**: The maximum number of PDB templates to use in modeling, `PDB_TO_USE`. This helps in narrowing down the best templates for accurate modeling.
- **Residues Cutoff**: A cutoff for residues, `MODEL_CUTOFF`, where non-covered regions in the alignment that are this length or longer will be replaced by chain breaks.
- **Input File**: Path to the input file containing mutation information, accessed via `args.input_file`.
- **Output Directory**: The directory where output files will be stored, `args.out_path`.
- **Databases Directory**: Location of the databases downloaded from Uniprot, including isoforms, specified in `DB_LOCATION`.
- **COBALT**: Path to the COBALT executable, `COBALT_HOME`, used for sequence alignment.
- **HMMER**: Location of the HMMER software executable, `HMMER_HOME`, essential for sequence analysis.
- **Weight Structural Score**: The weight assigned to the structural score in the modeling process, `W_STRUCT`. For details on how weights are utilized, refer to [Oliva et al.](https://doi.org/10.3389/fchem.2022.1059593).
- **Weight Mutation Score**: Specifies the weight for the mutation score, `W_MUT`, as described in the referenced publication by Oliva et al.
- **Modeller Executable**: Path to the MODELLER executable, `MODELLER_EXEC`, for homology or comparative modeling of protein three-dimensional structures.



## Class Methods

### `__init__(self, *args, **kwargs)`

Initializes the parser with any arguments and keyword arguments acceptable by `argparse.ArgumentParser`, and sets up the expected arguments.

### `add_arguments(self) -> None`

Defines and adds expected command line arguments to the parser. These include paths to input and output files, paths to external tools like COBALT and HMMER, and computational parameters such as resolution and sequence identity.

### `check_arguments(self, args: argparse.Namespace) -> None`

Checks the provided command line arguments to ensure that all necessary information is available. Raises an error if required arguments are missing or if there's a conflict between the parameters file and manual input.

### `parse_input(self, mutation_file_path: str) -> List`

Parses a given file containing a list of mutations and genes. Returns a structured list containing this information for further processing.

### `get_parameters(self, parameters_path: str = "parameters.dat") -> Dict`

Loads and returns parameters from a specified parameters file. Defaults to 'parameters.dat' if no file path is provided. This method structures the parameters into a dictionary for easy access.

### `merge_parameters(self, args: argparse.Namespace) -> Dict`

Merges parameters provided through the command line with those loaded from the parameters file. This ensures that command line arguments take precedence over file-based parameters.

### `print_parameters(self, args: argparse.Namespace, parameters: Dict) -> None`

Prints the final set of parameters that will be used for execution. This includes both command line arguments and parameters loaded from the file, highlighting the configuration for the user.

### `load_input(self, argv) -> tuple`

The main method for loading and validating user inputs. It parses command line arguments, validates them, loads additional parameters from a file if specified, and merges them as necessary. Returns a tuple containing the structured arguments and parameters.


