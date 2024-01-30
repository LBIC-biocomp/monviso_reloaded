import sys
from monviso_reloaded.parser import Parser

def main(argv=None) -> None:  # pragma no cover
    """
    Main function

    :param argv: argv

    :return: None
    """
    # arguments and parameters
    parser = Parser()
    args,parameters=parser.load_input(argv)


def init() -> None:
    if __name__ == "__main__":
        sys.exit(main())
init()
