from typing import Any, Dict, List, Union

class Gene:
    def __init__(self,gene_mutation_block: List[List]):
        self.name=gene_mutation_block[0]
        self.mutations=gene_mutation_block[1:]

class Isoform:
    def __init__(self):
        pass