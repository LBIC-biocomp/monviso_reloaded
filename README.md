
# Monviso reloaded
MoNvIso is a comprehensive software tool designed for the analysis and modeling of protein isoforms. It automates the process of identifying canonical and additional isoforms, assessing their modeling propensity, mapping mutations accurately, and building structural models of proteins. By leveraging data from the Uniprot database, MoNvIso facilitates a deeper understanding of protein function and variation.

[Read the documentation here.](https://alisamalb.github.io/monviso_reloaded/)

## Quick start

 - Download uniprot_sprot.fasta and uniprot_sprot_varsplic.fasta from UniProt.

 - Write a mutation.txt file with the gene name and the mutations you are interested in (e.g.,):
```
GRIN1
R       844     C
Ala     349     Thr
Pro     578     Arg
Ser     688     Tyr
Tyr     647     Ser

GRIN2B
E413G
C436R
M1342R
L1424F
PRO1439ALA
```

 - Write a parameters.dat with the calculation parameters:

```
DB_LOCATION=./
COBALT_HOME=./ncbi-cobalt-3.0.0/bin/
HMMER_HOME=/usr/local/bin/  
MODELLER_EXEC=mod10.5
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


```bash
$ monviso_reloaded -i mutations.txt -o out/ -p parameters.dat
```

## Development

Read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
