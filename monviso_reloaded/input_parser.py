import argparse
from pathlib import Path
from typing import Dict, List


class InputParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_arguments()

    def add_arguments(self) -> None:
        """Add a list of expected arguments to the Parser instance.

        :return: None
        """
        super().add_argument(
            "-i",
            "--input_file",
            help="Path to the file genes and mutations",
            type=str,
            required=True,
        )

        super().add_argument(
            "-o",
            "--out_path",
            help="Path to the output folder",
            type=str,
            required=True,
        )

        super().add_argument(
            "-pf",
            "--par_file",
            help="Input the path to the parameters file",
            type=str,
            required=False,
        )
        
        super().add_argument(
            "-pst",
            "--pesto_home",
            help="path to the directory containing the\
                PeSto repo downloaded from GitHub",
            type=str,
            required=False,
        )
        
        super().add_argument(
            "-db",
            "--db_home",
            help="path to the DBs containing the\
                canonical and isoform sequences",
            type=str,
            required=False,
        )

        super().add_argument(
            "-cobalt",
            "--cobalt_home",
            help="path to COBALT bin folder",
            type=str,
            required=False,
        )
        
        super().add_argument(
            "-ms",
            "--msms_home",
            help="Path to the directory containing the executables of MSMS",
            type=str,
            required=False,
        )

        super().add_argument(
            "-hmmer",
            "--hmmer_home",
            help="path to HMMER bin folder",
            type=str,
            required=False,
        )

        parameters_group = super().add_argument_group(
            "input", "manually provides the inputs"
        )
        parameters_group.add_argument(
            "-res",
            "--resolution",
            help="worst accepted resolution",
            type=float,
            default=4.50,
            required=False,
        )
        parameters_group.add_argument(
            "-seqid",
            "--sequence_identity",
            help="minimum sequence identity accepted",
            type=float,
            default=25,
            required=False,
        )
        parameters_group.add_argument(
            "-max_pdb",
            "--max_pdb_templates",
            help="maximum number of PDBs to use as templates",
            type=int,
            default=10,
            required=False,
        )

        parameters_group.add_argument(
            "-wt_mod",
            "--max_model_wt",
            help="maximum number of wild-type models to build",
            type=int,
            default=1,
            required=False,
        )

        parameters_group.add_argument(
            "-mut_mod",
            "--max_model_mut",
            help="maximum number of mutant models to build",
            type=int,
            default=1,
            required=False,
        )

        parameters_group.add_argument(
            "-rc",
            "--res_cutoff",
            help="maximum number of not crystalised residues\
                accepted before cutting a model",
            type=float,
            default=5.0,
            required=False,
        )

        parameters_group.add_argument(
            "-ws",
            "--w_struct",
            help="weight of the structural function",
            type=float,
            default=10.0,
            required=False,
        )

        parameters_group.add_argument(
            "-wm",
            "--w_mut",
            help="weight of the mutation function",
            type=float,
            default=10.0,
            required=False,
        )

        parameters_group.add_argument(
            "-mod",
            "--modeller_exec",
            help="path to the modeller executable",
            type=str,
            default="mod10.4",
            required=False,
        )

    def check_arguments(self, args: argparse.Namespace) -> None:
        """Verify that the necessary parameters have been provided

        :param args: passed arguments
        """
        if not args.par_file and (not args.db_home or not args.cobalt_home
                                  or not args.hmmer_home
                                  or not args.pesto_home):
            raise TypeError(
                "Specify a parameters file or insert the paths\
                    to the DBs, COBALT, PESTO and HMMER manually"
            )
        elif args.par_file and (
            args.db_home or args.cobalt_home or args.hmmer_home
        ):
            raise TypeError(
                "Either specify a parameters file or\
                    insert parameters manually"
            )

    def parse_input(self, mutation_file_path: argparse.Namespace) -> List:
        """Parse the list of mutations and genes from the mutation_list file.

        :param mutation_list: path to the file containing the list
        of mutations and genes
        :return: The list of gene and mutations
        """
        with Path(mutation_file_path).open() as my_file:
            content = my_file.read()
            
        
        #Remove useless new lines
        while content[-1]=="\n":
            content=content[:-1]
            
        while "\n\n\n" in content:
            content=content.replace('\n\n\n','\n\n')
        
        blocks = [block.splitlines() for block in content.split("\n\n")]
        return blocks
    
    def parse_sequences(self,sequence_file_path: argparse.Namespace) -> List:
        """Parses a sequence file to extract gene names, sequence names,
        mutation lists, and protein sequences.

        Args:
            sequence_file_path (argparse.Namespace): The path to the file containing
                                                    gene and sequence names, followed by
                                                    mutation list and on protein sequence
                                                    on a new line.

        Returns:
            List: An alternating list of [gene name, sequence name, mutation list] followed by
                                        the string of the sequence.
        """
        parsed_sequences=[]
        with Path(sequence_file_path).open() as my_file:
            content = my_file.read().split("\n")
        
        for i, line in enumerate(content):
            
            if line=="\n":
                pass
            elif line.startswith(">"):
                splitline=line[1:].split(':')
                if len(splitline)<1:
                    raise ValueError('While reading the sequence metadata, not enough information was found.'
                                     '\nThe first line must start with ">" and include at least a protein name'
                                     ' and a sequence name separated by ";".'
                                     f"The line that led to the error is: {line}")
                if len(splitline)>3:
                    raise ValueError('Too many columns (:) in the file.'
                                     f"The line that led to the error is: {line}")
                if len(splitline)==3:
                    splitline[2]=[mutation.replace(" ","") for mutation in splitline[2].split(",")]
                    
                parsed_sequences.append([splitline])            
            else:
                parsed_sequences[-1].append(line)
        
        if min([len(chunk) for chunk in parsed_sequences])<2:
            raise ValueError("Not all lines starting with > have an associated sequence.")      
        else:
           merged_parsed_sequences=[]
           for chunk in parsed_sequences:
               merged_parsed_sequences.append(chunk[0])
               merged_parsed_sequences.append("".join(chunk[1:]).replace("\n",""))
        
        return merged_parsed_sequences

    def get_parameters(
        self,
        parameters_path: argparse.Namespace = "parameters.dat",
    ) -> Dict:
        """Collect the parameters from the parameters file if provided.

        :param parameters_path: Path to the parameters file,
        defaults to "parameters.dat".

        :return: Dict of keywords with the associated parameters
        """
        keywords = {
            "RESOLUTION": None,
            "SEQID": None,
            "PDB_TO_USE": None,
            "PESTO_HOME": None,
            "INPUT_FILE": None,
            "OUTPUT_PATH": None,
            "DB_LOCATION": None,
            "HMMER_HOME": None,
            "COBALT_HOME": None,
            "MSMS_HOME": None,
            "MODEL_CUTOFF": None,
            "NUM_OF_MOD_WT": None,
            "NUM_OF_MOD_MUT": None,
            "W_STRUCT": None,
            "W_MUT": None,
            "MODELLER_EXEC": None,
        }
        keys = list(keywords)
        if not Path(parameters_path).exists():

            error_message = f"Parameters file not found in {parameters_path},\
                please check the path provided."
            raise TypeError(error_message)

        with Path(parameters_path).open() as my_file:
            lines = my_file.readlines()

        for key in keys:
            for line in lines:
                if key in line:
                    value = line[line.find("=") + 1 :].strip()
                    keywords[key] = value
        return keywords

    def convertArgsToParameters(self, args: argparse.Namespace) -> Dict:
        """Make the parameters dictionary if they are passed as arguments

        :param args: passed arguments
        """
        
        return {
            "RESOLUTION": args.resolution,
            "SEQID": args.sequence_identity,
            "PDB_TO_USE": args.max_pdb_templates,
            "PESTO_HOME":args.pesto_home,
            "OUTPUT_PATH":args.out_path,
            "INPUT_FILE": args.input_file,
            "DB_LOCATION": args.db_home,
            "HMMER_HOME": args.hmmer_home,
            "COBALT_HOME": args.cobalt_home,
            "MSMS_HOME": args.msms_home,
            "MODEL_CUTOFF": args.res_cutoff,
            "NUM_OF_MOD_WT": args.max_model_wt,
            "NUM_OF_MOD_MUT": args.max_model_mut,
            "W_STRUCT": args.w_struct,
            "W_MUT": args.w_mut,
            "MODELLER_EXEC": args.modeller_exec,
        }
        
    def merge_parameters(self,args,parameters):
        """Merges args and parameters, with a preference for parameters.
        
        The settings expressed in the command line will be added to the parameters dictionary.
        """
                
        args_name=[args.resolution,
              args.sequence_identity,
              args.max_pdb_templates,
              args.pesto_home,
              args.input_file,
              args.out_path,
              args.db_home,
              args.hmmer_home,
              args.cobalt_home,
              args.msms_home,
              args.res_cutoff,
              args.max_model_wt,
              args.max_model_mut,
              args.w_struct,
              args.w_mut,
              args.modeller_exec
              ]
        
        dict_keys=['RESOLUTION',
                   'SEQID',
                   'PDB_TO_USE',
                   'PESTO_HOME',
                   'INPUT_FILE',
                   'OUTPUT_PATH',
                   'DB_LOCATION',
                   'HMMER_HOME',
                   'COBALT_HOME',
                   'MSMS_HOME',
                   'MODEL_CUTOFF',
                   'NUM_OF_MOD_WT',
                   'NUM_OF_MOD_MUT',
                   'W_STRUCT',
                   'W_MUT',
                   'MODELLER_EXEC']
        
        
        for i, key in enumerate(parameters):
            if parameters[key] is None:
                try:
                    parameters[key]=args_name[i]
                except:
                    raise RuntimeError(f"The argument for {dict_keys[i]} was\
                    not found in the parameter file,\
                    and there was a error in reading it from the\
                    command line .")
                    
    def print_parameters(
        self, args: argparse.Namespace, parameters: Dict
    ) -> None:
        """Print the parameters provided

        :param parameters: Provided parameters
        """
        param = (
            f"\nResolution: {parameters['RESOLUTION']}\n"
            f"SEQ ID: {parameters['SEQID']}\n"
            f"WT MODELS TO PREPARE: {parameters['NUM_OF_MOD_WT']}\n"
            f"MUTANTS MODEL TO PREPARE: {parameters['NUM_OF_MOD_MUT']}\n"
            f"MAX PDBS AS TEMPLATES: {parameters['PDB_TO_USE']}\n"
            f"RESIDUES CUTOFF: {parameters['MODEL_CUTOFF']}\n"
            f"INPUT FILE: {args.input_file}\n"
            f"OUTPUT DIRECTORY: {args.out_path}\n"
            f"DATABASES DIRECTORY: {parameters['DB_LOCATION']}\n"
            f"COBALT: {parameters['COBALT_HOME']}\n"
            f"PESTO: {parameters['PESTO_HOME']}\n"
            f"HMMER: {parameters['HMMER_HOME']}\n"
            f"WEIGHT STRUCTURAL SCORE: {parameters['W_STRUCT']}\n"
            f"WEIGHT MUTATION SCORE: {parameters['W_MUT']}\n"
            f"MODELLER EXECUTABLE: {parameters['MODELLER_EXEC']}"
        )

        print(param)

    def load_input(self, argv) -> tuple:
        """Load user input from the command line and parameters file.

        :param argv: command line arguments
        """
        args, unparsed = self.parse_known_args(argv)
        self.check_arguments(args)
        if args.par_file:
            parameters = self.get_parameters(args.par_file)
        else:
            parameters = self.convertArgsToParameters(args)
        
        self.merge_parameters(args,parameters)
        self.print_parameters(args, parameters)
        return (args, parameters)
