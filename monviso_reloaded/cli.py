import sys

from monviso_reloaded.base import IsoformRun
from monviso_reloaded.base import AnalysisRun

def isoform(argv):
    run = IsoformRun()
    run.load_input(argv[1:])
    run.load_mutation_list()
    run.create_genes()
    run.create_isoforms()
    run.run_blastp()
    run.run_cobalt()
    run.run_hmmsearch()
    run.load_templates()
    run.select_isoforms()
    run.start_modeller()
    run.write_report()

def analysis(argv):
    run = AnalysisRun()
    run.load_input(argv[1:])
    run.load_genes_from_mutation_list()
    run.pestoAnalysis()
    
def main(argv=None):  # pragma no cover
    """
    Main function

    :param argv: argv

    :return: None
    """
    # arguments and parameters
    if argv is None:
        argv = sys.argv[1:]
    
    if argv[0]=="isoform":
        isoform(argv)

    elif argv[0]=="analysis":
        analysis(argv)
    
    else:
        print("The command was not recognized. The available commands are:")
        print(" - monviso isoform")
        print(" - monviso analysis")

def init() -> None:
    if __name__ == "__main__":
        sys.exit(main())


init()
