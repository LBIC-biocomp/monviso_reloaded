import sys

from monviso_reloaded.base import HomologyModelRun
from monviso_reloaded.base import AnalysisRun
from monviso_reloaded.base import DockingRun


def isoform(argv):
    run = HomologyModelRun()
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

def sequence(argv):
    run = HomologyModelRun()
    run.load_input(argv[1:])
    run.load_sequence_list()
    run.create_genes(sequenceRun=True)
    run.create_isoforms(sequenceRun=True)
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
    run.analysis()

def dock(argv):
    run = DockingRun()
    run.load_input(argv[1:])
    run.load_genes_from_list()
    run.analysis()

def unrecognized_command():
    print("The command was not recognized. The available commands are:")
    print(" - monviso isoform")
    print(" - monviso sequence")
    print(" - monviso analysis")
    print(" - monviso dock")
    print("append \"--help\" at the end to receive help on the tool.")
    
def main(argv=None):  # pragma no cover
    """
    Main function

    :param argv: argv

    :return: None
    """
    # arguments and parameters
    if argv is None:
        argv = sys.argv[1:]
    
    if len(argv)>0:
        if argv[0]=="isoform":
            isoform(argv)

        elif argv[0]=="analysis":
            analysis(argv)
            
        elif argv[0]=="sequence":
            sequence(argv)
            
        elif argv[0]=="dock":
            dock(argv)
        else:
            unrecognized_command()
    else:
        unrecognized_command()

def init() -> None:
    if __name__ == "__main__":
        sys.exit(main())


init()
