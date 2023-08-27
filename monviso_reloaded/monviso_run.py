from utils.utils import parse_input
from utils.utils import get_parameters
from utils.utils import make_gene_directories


mut_file = "/mnt/d/work/monviso/mutations.txt"
param_file = "/mnt/d/work/monviso/parameters.dat"
#output_dir = "/mnt/d/work/monviso/parameters.dat"


blocks, protein_list = parse_input(mut_file)

parameters = get_parameters(param_file)   # processes the parameters
print(parameters)

make_gene_directories(blocks, parameters["OUT_DIR"])