from typing import Any, Dict, List, Union
from .parser import Parser
from .file_handler import FileHandler
from .gene import Gene

class Run:
    def __init__(self):
        self.args=[]
        self.parameters=[]
        self.mutation_list=[]
        self.parser=Parser()
        self.genes=[]
    
    def load_input(self,argv) -> None:
        """Load user input from the command line and parameters file
         and save them as attributes.

        :param argv: command line arguments
        """
        self.args,self.parameters=self.parser.load_input(argv)
        
    def load_mutation_list(self) -> None:
        """Parse the list of mutations and genes from the mutation_list file
        and save it as attribute.

        :param mutation_list: path to the file containing the list of mutations and genes
        """
        self.mutation_list=self.parser.parse_input(self.args.input_file)
        
    def create_genes(self) -> None:
        """Take the list of mutations of the genes saved in self.mutation_list, and
        for each gene, create a Gene instance. All Gene instances are saved in the
        self.gene list.

        :param mutation_list: A list of the blocks extracted from the input file.
        At index 0, the list contains the name of the gene.
        """
        for i,gene_mutation_block in enumerate(self.mutation_list):
            self.genes.append(Gene(gene_mutation_block))
