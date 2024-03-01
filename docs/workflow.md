# Workflow

## For the user
- The user provides the parameters (software paths, model settings, ...)
- The user provides a list of genes and mutations
- MoNvIso searches UniProt's databases for the canonical sequence and other isoforms
- In the output path, a directory is created for every gene
- In every gene's directory, a sub-directory is created for every isoform
- The PDB databank is queried for every sequence and returns a list of similar sequences
- HMMER is used to combine the results and start a new search for templates
- The PDB files of the templates are downloaded, and the template chain is extracted
- The sequence of the extracted chains is aligned with the sequence to model
- Chain breaks are introduced in place of large gaps
- The number of mutations that are mappable and the value of sequence identity are measured
- Isoforms are ranked based on their score (ability to model the mutated protein)
- The aligned sequence is passed to Modeller for modelling

## Under the hood
The MoNvIso run is coordinated by an instance of the Run class.

- The `Run.load_mutation_list()` method passes the mutation file to the Input parser, which returns the mutations in a list type
- The `Run.load_input()` method takes the input parameters and passes them to the Input_parser file, which returns them processed as a dictionary
- The `Run.create_genes()` method instances the Gene object, one object per gene to model. These genes are stored in a list as a local attribute of the Run.
- The `Run.create_isoform()` method loads the UniProt databases via a Database_Parser object. The database returns one sequence per isoform per gene. Each is used to instance an Isoform object. Every isoform object is stored in a list as an attribute of the corresponding gene.
- At this point, the Gene object and Isoform objects take care of creating their corresponding directories and .fasta files with the alignments
- The Run object coordinates the creation of alignments with the `Run.run_blastp()`, `Run.run_cobalt()`, `Run.run_hmmsearch()` methods. Each isoform of each gene runs a blastp search, saves the output, creates a .hmm with HMMER, and uses it for a search with hmmsearch. The output contains the PDB IDs and chains of the best templates
- `Run.load_templates()`: For every proposed template, a Template instance is created, which takes care of downloading and extracting the chain and fasta sequence from the PDB file. The operations on the PDB file are managed by the PDB_manager class. All the sequences of the templates are aligned with Cobalt.
- `Run.select_isoform()`: For every isoform, the mutation score (How many mutations can be mapped on the sequence), structural score (sequence identity), and the Selection score (the first two combined) are calculated. The isoforms are ranked, and the ones with higher scores are preferred to model the mutation.
- `Run.start_modeller()`: The isoforms that can be modeled will create the Modeller input file and start a Modeller run.
