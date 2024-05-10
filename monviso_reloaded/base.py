from .database_parser import DatabaseParser
from .gene import Gene
from .input_parser import InputParser
from .analyzer import Analyzer

class IsoformRun:
    def __init__(self):
        self.parameters = []
        self.mutation_list = []
        self.input_parser = InputParser()
        self.genes = []

    def load_input(self, argv) -> None:
        """Load user input from the command line and parameters file
         and save them as attributes.

        :param argv: command line arguments
        """
        _, self.parameters = self.input_parser.load_input(argv)

    def load_mutation_list(self) -> None:
        """Parse the list of mutations and genes from the mutation_list file
        and save it as attribute.
        """
        self.mutation_list = self.input_parser.parse_input(
            self.parameters["INPUT_FILE"]
        )

    def create_genes(self) -> None:
        """Take the list of mutations of the genes saved
        in self.mutation_list, and for each gene, create a Gene instance.
        All Gene instances are saved in the
        self.gene list.
        """
        for i, gene_mutation_block in enumerate(self.mutation_list):
            self.genes.append(Gene(gene_mutation_block, self.parameters["OUTPUT_PATH"]))

    def create_isoforms(self) -> None:
        """Loads isoforms for each gene in the 'genes' attribute from
        the Uniprot database.
        """
        with DatabaseParser(self.parameters["DB_LOCATION"]) as db_parser:
            for gene in self.genes:
                gene.load_isoforms(db_parser)

    def run_blastp(self) -> None:
        """Start a blastp query for every loaded isoform. The proper
        function is a method of the class Isoform. This is used to
        split the whole in smaller steps.
        """
        for gene in self.genes:
            for isoform in gene.isoforms:
                isoform.blastp_search()

    def run_cobalt(self) -> None:
        """Start a cobalt run for every loaded isoform. The proper
        function is a method of the class Isoform. This is used to
        split the whole in smaller steps.
        """
        for gene in self.genes:
            for isoform in gene.isoforms:
                isoform.create_MSA(cobalt_home=self.parameters["COBALT_HOME"])

    def run_hmmsearch(self) -> None:
        """Start a hmmsearch for every loaded isoform. The proper
        function is a method of the class Isoform. This is used to
        split the whole in smaller steps.
        """
        for gene in self.genes:
            for isoform in gene.isoforms:
                isoform.HMMsearch(hmmer_home=self.parameters["HMMER_HOME"])

    def load_templates(self) -> None:
        for gene in self.genes:
            for isoform in gene.isoforms:
                isoform.load_templates(
                    int(self.parameters["PDB_TO_USE"]),
                    float(self.parameters["RESOLUTION"]),
                    self.parameters["COBALT_HOME"],
                )

    def select_isoforms(self) -> None:
        for gene in self.genes:
            gene.select_isoforms(
                float(self.parameters["W_STRUCT"]),
                float(self.parameters["W_MUT"]),
                float(self.parameters["SEQID"]),
                int(self.parameters["MODEL_CUTOFF"]),
                cobalt_home=self.parameters["COBALT_HOME"],
            )

    def start_modeller(self) -> None:
        """Run modeller only for the isoforms to model.
        The method write_modeller() of the isoform accepts the mutation
        as argument.
        """
        for gene in self.genes:
            gene.report_on_selected_isoforms()
            for isoform, mutation in gene.isoforms_to_model:
                isoform.run_modeller(
                    mutation, self.parameters["MODELLER_EXEC"],
                    int(self.parameters["MODEL_CUTOFF"]),
                    int(self.parameters["NUM_OF_MOD_WT"]),
                    int(self.parameters["NUM_OF_MOD_MUT"]),
                    self.parameters["COBALT_HOME"])

            

    def write_report(self):
        for gene in self.genes:
            gene.write_report()
            

class AnalysisRun:
    def __init__(self):
        self.parameters = []
        self.genes = []
        self.input_parser = InputParser()
        self.genes = []
        
    def load_input(self, argv) -> None:
        """Load user input from the command line and parameters file
         and save them as attributes.

        :param argv: command line arguments
        """
        _, self.parameters = self.input_parser.load_input(argv)
        self.pesto_home=self.parameters["PESTO_HOME"]
        self.output_path=self.parameters["OUTPUT_PATH"]

    def load_genes_from_mutation_list(self) -> None:
        """Parse the list of mutations and genes from the mutation_list file
        and save the genes as attribute.
        """
        self.genes = [m[0] for m in self.input_parser.parse_input(
            self.parameters["INPUT_FILE"])]
        
    def analysis(self) -> None:
        self.analyzer=Analyzer(self.pesto_home,self.output_path,
                                    self.genes)
        
