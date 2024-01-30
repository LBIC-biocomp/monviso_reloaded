from .parser import Parser

class Run:
    def __init__(self):
        self.args=[]
        self.parameters=[]
        self.mutation_list=[]
        self.parser=Parser()
    
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
        blocks,protein_list=self.parser.parse_input(self.args.input_file)
        print(blocks)
        print(protein_list)
        


class Gene:
    def __init__(self):
        pass

class Isoform:
    def __init__(self):
        pass