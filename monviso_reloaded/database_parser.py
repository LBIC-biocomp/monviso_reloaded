from pathlib import Path


class DatabaseParser:
    def __init__(self, db_location: str):
        self.db_location = db_location
        self.load_database()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        pass

    def parse_database(self, database_path: Path) -> list[list[str]]:
        """
        Splits the sequences from a database file into a list, where each
        list element corresponds to a single sequence.
        Each sequence is represented as a list of strings, with each string
        being a line from the database that contains
        a part of the sequence. This allows for multi-line sequences in the
        database to be captured fully in each element.

        :param database_path: The Path pointing to the Uniprot database
        :return: A list of gene sequences, each split in multiple lines
        """
        with open(database_path, "r") as f:
            content = f.read()
        return [
            block.splitlines()
            for block in content.split(">")
            if "sp|" in block
        ]

    def load_database(self):
        """
        Load Uniprot databases as local attributes, using the
        parse_database() method
        """
        canonical_db_path = Path(self.db_location, "uniprot_sprot.fasta")
        print("Loading canonical database...")
        self.canonical_db = self.parse_database(canonical_db_path)
        print(f"Loaded {len(self.canonical_db)} elements.")

        print("Loading split variants database...")
        isoforms_db_path = Path(
            self.db_location, "uniprot_sprot_varsplic.fasta"
        )
        self.isoforms_db = self.parse_database(isoforms_db_path)
        print(f"Loaded {len(self.isoforms_db)} elements.")

    def get_canonical_isoforms(self, gene_name: str) -> list[list[str]]:
        """Use the gene name to retrieve the corresponding sequences
        from the canonical database.

        :param gene_name: Name of the gene
        :return: List of multi-line gene sequences
        """
        gene_name = f"GN={str(gene_name.upper())}"
        species = "OS=Homo sapiens"
        output = [seq for seq in self.canonical_db if gene_name in seq[0] and species in seq[0]]
        if len(output) == 0:
            print(
                f"No sequences found for gene {gene_name} in\
the canonical isoforms database."
            )
        return output

    def get_noncanonical_isoforms(self, gene_name: str) -> list[list[str]]:
        """Use the gene name to retrieve the corresponding sequences from
        the split variants database.

        :param gene_name: Name of the gene
        :return: List of multi-line gene sequences
        """
        gene_name = f"GN={str(gene_name.upper())}"
        species = "OS=Homo sapiens"
        output = [seq for seq in self.isoforms_db if gene_name in seq[0] and species in seq[0]]
        if len(output) == 0:
            print(
                f"No sequences found for gene {gene_name} in \
the split variants isoforms database."
            )
        return output
