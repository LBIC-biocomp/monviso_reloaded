import argparse
import sys
from pathlib import Path
import os

from monviso_reloaded.utils.utils import add_arguments
from monviso_reloaded.utils.utils import check_arguments
from monviso_reloaded.utils.utils import get_parameters
from monviso_reloaded.utils.utils import print_parameters
from monviso_reloaded.utils.utils import parse_input
from monviso_reloaded.utils.utils import make_gene_directories
from monviso_reloaded.run.preparation import get_isoforms_from_db
from monviso_reloaded.run.preparation import build_master_isoform_file


#from monviso_reloaded.utils.utils import merge_parameters


#from monviso_reloaded.run.file1 import run_hmm



# mut_file = "/mnt/d/work/monviso/mutations.txt"
# param_file = "/mnt/d/work/monviso/parameters.dat"
#output_dir = "/mnt/d/work/monviso/parameters.dat"

def main(argv=None) -> None: #pragma no cover
    """
    Main function

    :param argv: argv

    :return: None
    """
    # arguments and parameters
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    args, unparsed = parser.parse_known_args(argv)
    check_arguments(args)
    if args.par_file:
        parameters = get_parameters(args.par_file)
    else:
        parameters = merge_parameters(args)  
              
    print_parameters(args, parameters)

    blocks, protein_list = parse_input(args.input_file)
    print(protein_list)
    make_gene_directories(blocks, args.out_path)
    get_isoforms_from_db(protein_list, args.out_path, 
                             parameters["DB_LOCATION"])

    for gene in protein_list:
        gene_path = Path(args.out_path, gene)
        os.chdir(gene_path) 
        # Changing directory (i.e., previous line) not compatible with pathlib
	# gene_path cannot be used as argument for build_baster_isoform_file
        build_master_isoform_file(Path(os.getcwd()))
        check_output = run_hmm(gene_path, gene, parameters["RESOLUTION"],
                            parameters["SEQID"], parameters["HMMER_HOME"],
                            parameters["COBALT_HOME"], parameters["PDB_TO_USE"])



def init() -> None:
    if __name__ == "__main__":
        sys.exit(main())


init()
