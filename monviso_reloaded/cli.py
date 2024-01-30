import sys
from monviso_reloaded.base import Run

def main(argv=None) -> None:  # pragma no cover
    """
    Main function

    :param argv: argv

    :return: None
    """
    # arguments and parameters
    run= Run()
    run.load_input(argv)
    run.load_mutation_list()


def init() -> None:
    if __name__ == "__main__":
        sys.exit(main())
init()
