# Cobalt Wrapper Class

## Overview

The `Cobalt` class serves as a Python wrapper for COBALT (Constraint-based Multiple Alignment Tool), a software tool used for multiple sequence alignment (MSA). This class enables the integration of COBALT's alignment capabilities within Python applications, facilitating the automated processing of sequence alignments. 

## Methods

### `__init__(self)`

Initializes a new instance of the `Cobalt` class. This method prepares the environment for running COBALT alignments but does not require any parameters.

### `run(self, hits_file: Union[str, Path], aligned_file: Union[str, Path], cobalt_home: str) -> bool`

Executes the COBALT alignment process on a specified set of sequences.

- **Parameters:**
  - `hits_file`: The path to the input file containing sequence data to be aligned, typically the output from a BLASTP search. This can be a `str` or `Path` object.
  - `aligned_file`: The path to the output file where the aligned sequences will be saved, in multi-FASTA format. This can be a `str` or `Path` object.
  - `cobalt_home`: The file system path to the COBALT installation directory, which contains the COBALT executable.

- **Returns:**
  - `bool`: Returns `True` if the COBALT alignment completes successfully, otherwise returns `False`.
