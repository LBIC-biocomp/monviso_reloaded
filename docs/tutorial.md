# Quick Start Tutorial for MonViso Reloaded

This tutorial will walk you through the process of preparing input files and running MonViso Reloaded for mutation analysis.

## Step 1: Download the UniProt databases

Download and decompress the following files:

- [Database of isoforms](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz)
- [Database of split variances](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz)

Remember the directory containing the .fasta files. It will be one of the parameters required to run the software.

## Step 2: Prepare Input Files

You will need to create two files: `mutations.txt` for listing the mutations and `parameters.dat` for specifying various parameters and paths used by MonViso Reloaded.

### Writing `mutations.txt`

Create a text file named `mutations.txt` and fill it with your mutations of interest, formatted as shown below:

```text
GRIN1
R 844 C
Ala 349 Thr
Pro 578 Arg
Ser 688 Tyr
Tyr 647 Ser

GRIN2B
E413G
C436R
M1342R
L1424F
PRO1439ALA
```


Each section starts with the gene name, followed by lines specifying the original amino acid, its position, and the mutated amino acid.

### Writing `parameters.dat`

Create a text file named `parameters.dat` and specify the configuration for your analysis as follows:

```
DB_LOCATION= [path to database location]
COBALT_HOME=[path to cobalt binaries]
HMMER_HOME=[path to HMMER binaries]
MODELLER_EXEC=[command to run modeller e.g., mod10.5]
RESOLUTION=4.50
SEQID=25
HMM_TO_IMPORT=100
MODEL_CUTOFF=5
PDB_TO_USE=10
NUM_OF_MOD_WT=1
NUM_OF_MOD_MUT=1
W_STRUCT=10
W_MUT=10
```

## Step 3:   Run the command

```bash
monviso_reloaded -i mutations.txt -o out -p parameters.dat
```