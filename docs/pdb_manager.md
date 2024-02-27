# PDB Manager

The `PDBManager` class is a tool designed for managing Protein Data Bank (PDB) files in Monviso. It leverages Biopython's capabilities to download, filter, and extract valuable information from PDB files. 

## Features

- **Download PDB Files**: Download PDB files from the Protein Data Bank.
- **Extract Clean Chain**: Extracts a specific chain from a PDB file, filtering out non-standard residues and atoms, and checks for resolution cutoffs.
- **Extract FASTA Sequence**: Extracts the FASTA sequence from a PDB file.

## `PDBManager`

Manages the downloading, cleaning, and processing of PDB files.

### Methods

- `downloadPDB(pdb, out_path)`: Downloads a PDB file and saves it locally.
- `extract_clean_chain(input_pdb_path, output_pdb_path, chain_letter, resolution_cutoff)`: Extracts a specific chain from a PDB file, ensuring it meets resolution quality standards.
- `extract_fasta(pdb_name, pdb_path, output_fasta_path)`: Extracts the FASTA sequence from a PDB file.


## Accessory class

### `ChainSelection(Select)`

A subclass of Biopython's `Select` class for filtering specific chains and residues in a PDB file.

#### Methods

- `accept_model(model)`: Accepts only the first model in a PDB file.
- `accept_residue(residue)`: Filters residues to include only standard amino acids.
- `accept_chain(chain)`: Filters for a specific chain.
- `accept_atom(atom)`: Filters for standard atoms.
