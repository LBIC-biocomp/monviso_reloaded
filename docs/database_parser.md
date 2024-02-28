# Database Parser Class

## Overview

The `DatabaseParser` class is designed for processing Uniprot databases to retrieve sequences for specific genes. It facilitates the parsing of canonical and noncanonical (split variants) isoform databases downloaded from Uniprot. 

## Methods

### `__init__(self, db_location: str)`

Initializes the `DatabaseParser` with the location of the Uniprot databases. It automatically loads the databases upon initialization.

- **Parameters:**
  - `db_location`: A string specifying the filesystem path to the directory containing the Uniprot database files.

### `parse_database(database_path: Path) -> list[list[str]]`

Parses a given database file, organizing its content into a structured list format where each element represents a single sequence.

- **Parameters:**
  - `database_path`: A `Path` object pointing to the Uniprot database file to be parsed.
- **Returns:**
  - A list of gene sequences, each represented as a list of strings corresponding to lines from the database.

### `load_database()`

Loads the canonical and noncanonical isoform databases from Uniprot, making them available as local attributes for query operations.

### `get_canonical_isoforms(gene_name: str) -> list[list[str]]`

Retrieves sequences corresponding to the specified gene from the canonical database.

- **Parameters:**
  - `gene_name`: The name of the gene for which sequences are being retrieved.
- **Returns:**
  - A list of sequences for the specified gene, with each sequence split across multiple lines.

### `get_noncanonical_isoforms(gene_name: str) -> list[list[str]]`

Retrieves sequences corresponding to the specified gene from the noncanonical (split variants) database.

- **Parameters:**
  - `gene_name`: The name of the gene for which sequences are being retrieved.
- **Returns:**
  - A list of sequences for the specified gene, with each sequence split across multiple lines.
