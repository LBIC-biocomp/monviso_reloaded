import sys

from monviso_reloaded.base import Run


def main(argv=None):  # pragma no cover
    """
    Main function

    :param argv: argv

    :return: None
    """
    # arguments and parameters
    run = Run()
    run.load_input(argv)
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


def init() -> None:
    if __name__ == "__main__":
        sys.exit(main())


init()
