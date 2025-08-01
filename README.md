![Monviso](https://github.com/LBIC-biocomp/monviso_reloaded/blob/main/logo.png?raw=true)

MoNvIso is a comprehensive software tool designed for the analysis and modeling of protein isoforms and protein-protein complexes. It automates the process of identifying canonical and additional isoforms, assessing their modeling propensity, mapping mutations accurately, and building structural models of proteins. The starting point can be either gene names or sequences, making it versatile for various workflows.

MoNvIso also automates SASA (solvent-accessible surface area) analysis, residue placement, and protein-protein interaction (PPI) identification. Furthermore, it streamlines protein-protein docking by combining the strengths of MEGADOCK, HDOCKLite, and HADDOCK 3, offering a seamless approach to complex structural predictions.



[Read the documentation here.](https://lbic-biocomp.github.io/monviso_reloaded/)

---

## Installation

1. Obtain a Modeller [license](https://salilab.org/modeller/registration.html).
2. Obtain an archive of [HDockLite](https://huanglab.phys.hust.edu.cn/software/hdocklite/).

3. (a) Run `installer.sh` on Linux:
   ```bash
   ./installer.sh
   ```
   
   (b) Alternatively, use the `Containerfile` with Docker:
   ```bash
   docker build --build-arg MODELLER_LICENSE=INSERT_LICENSE --build-arg HDOCKLITE_URL=PATH_TO_HDOCKlite.tar.gz -t monviso -f ./Containerfile ./

   ```

   (c) Install [msms](https://ccsb.scripps.edu/msms/), [HADDOCK3.0](https://github.com/haddocking/haddock3) with patched CNS1.3, [MEGADOCK](https://github.com/akiyamalab/MEGADOCK),       [COBALT](ftp://ftp.ncbi.nlm.nih.gov/pub/cobalt/executables/LATEST/), [HMMER](http://eddylab.org/software/hmmer/hmmer.tar.gz), [PeSTo](https://github.com/LBM-EPFL/PeSTo) manually.      Then clone this repository and
   ```bash
   cd monviso_reloaded/
   python -m pip install --upgrade pip setuptools wheel
   pip install -e .
   ```
   or
   ```bash
   pip install monviso
   ```



## Quick Setup

#### Modelling mutations from gene names
You will need:

- A file containing paths and parameters (parameters.dat),
- Uniprot database files and PDB databank sequences,
- A file listing gene names and mutations.

The parameter file will look like this:
```text
DB_LOCATION=/work/run
MEGADOCK_HOME=/work/Monviso/MEGADOCK
HDOCKLITE_HOME=/work/Monviso/HDOCKlite-v1.1
HADDOCK_HOME=/work/Monviso/haddock3
HADDOCK_SELECTION=sasa
HADDOCK_CUTOFF=0.9
MSMS_HOME=/work/Monviso/msms
COBALT_HOME=/work/Monviso/cobalt/bin
HMMER_HOME=/work/Monviso/hmmer/bin
PESTO_HOME=/work/Monviso/PeSTo
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

[`uniprot_sprot.fasta`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz ) and [`uniprot_sprot_varsplic.fasta`](https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz) must be downloaded from UniProt and placed in the DB_LOCATION path. The same for [`pdb_seqres.txt`](https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz), downloaded from RCSB.org and placed in the DB_LOCATION path. 

The mutation file lists the gene names and associated mutations. Example (mutations.txt):
```text
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

#### Modelling mutations from sequences
For modeling from sequences, you will need the UniProt databases, the parameters file, and a FASTA-format file with the sequences. Example (sequences.txt):
```text
>PROJ1:seq1:R3A,T7C
MKRISTTITTTITITTGNGAG

>PROJ2:seq2:
MPRVLLISLVTAALAAGQGDS
```


### Running Monviso

#### Isoform Modeling
Analyze isoforms and model structures from mutations:
```bash
monviso isoform -i mutations.txt -o out/ -pf parameters.dat
```

#### Modeling from Sequences
Model structures directly from sequences:
```bash
monviso sequence -i sequences.txt -o out/ -pf parameters.dat
```

#### SASA and PPI Interface Prediction
Perform SASA and residue depth analysis and PPI predictions with PeSTo:
```bash
monviso analysis -i mutations.txt -o out/ -pf parameters.dat
```

#### Protein-Protein Docking
Dock two modeled genes or sequences:
```bash
monviso dock -i dock.txt -o out/ -pf parameters.dat
```
Example `dock.txt`:
```text
GRIN1
GRIN2B

PROJ1
PROJ2
```

---


---

## Development

For contributing, read the [CONTRIBUTING.md](CONTRIBUTING.md) file.

---

