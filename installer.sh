#!/bin/bash

# Set up a local installation directory and variables
INSTALL_DIR="/work/Monviso"
MODELLER_LICENSE="MODELLER_LICENSE"
HDOCKLITE_URL="PATH_TO_HDOCKlite.tar.gz"

# Create the installation directory
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR

# Ensure necessary packages are installed (note: this might require sudo)
# If not possible, advise users to install them.
echo "Ensure wget, curl, git, build-essential are installed."


# Fail if the HDOCKLITE_URL is not provided
if [ -z "$HDOCKLITE_URL" ]; then
    echo "Error: HDOCKLITE_URL not provided"
    exit 1
fi


# Check if the URL is from 'huanglab' and download or copy accordingly
if echo "$HDOCKLITE_URL" | grep -q 'huanglab'; then
    echo "Downloading HDOCKlite from $HDOCKLITE_URL..."
    curl -o $INSTALL_DIR"/HDOCKlite.tar.gz" "$HDOCKLITE_URL"
else
    echo "Copying HDOCKlite from local path $HDOCKLITE_URL..."
    cp "$HDOCKLITE_URL" $INSTALL_DIR"/HDOCKlite.tar.gz"
fi

# Extract the downloaded or copied archive
echo "Extracting HDOCKlite..."
tar -xzf $INSTALL_DIR"/HDOCKlite.tar.gz" -C "$INSTALL_DIR"
mv HDOCK* HDOCKlite

# Clean up the archive
#rm $INSTALL_DIR"/HDOCKlite.tar.gz"

echo "HDOCKlite installation completed in $DEST_DIR."


# Download and install Miniconda locally
if ! command -v conda >/dev/null 2>&1; then
    echo "Downloading and installing Miniconda..."
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
    bash miniconda.sh -b -p $INSTALL_DIR/miniconda
    rm miniconda.sh
    export PATH="$INSTALL_DIR/miniconda/bin:$PATH"
else
    echo "Miniconda already installed."
fi

# Activate Conda and install Python environment
conda create -n monviso python=3.9 -y

# Install Modeller (requires MODELLER_LICENSE)
echo "Installing Modeller..."
/bin/bash -c "source activate monviso && conda install -c salilab modeller -y"

# Replace license in Modeller config (replace with your license)
MOD_PATH=$(which mod10.5)
if [ -z "$MOD_PATH" ]; then
    echo "mod10.5 not found. The installation of modeller did not work."
    exit 1
fi
BASE_PATH=$(dirname $(dirname "$MOD_PATH"))
CONFIG_PATH="$BASE_PATH/lib/modeller-10.5/modlib/modeller/config.py"
sed -i "s/XXXX/${MODELLER_LICENSE}/" $CONFIG_PATH

# Install Monviso using pip in the myenv environment
git clone https://github.com/LBIC-biocomp/monviso_reloaded
cd monviso_reloaded
/bin/bash -c "source activate monviso && pip install -e ."
cd ..

# Download and install HMMER locally
echo "Installing HMMER..."
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
tar xvzf hmmer.tar.gz
rm hmmer.tar.gz
mv hmmer-* source_hmmer
cd source_hmmer
./configure --prefix=$INSTALL_DIR/hmmer/
make
make install
cd ..
rm -rf source_hmmer

# Download and extract Cobalt
echo "Installing Cobalt..."
wget -r ftp://ftp.ncbi.nlm.nih.gov/pub/cobalt/executables/LATEST/*x64-linux.tar.gz
mv $INSTALL_DIR ftp.ncbi.nlm.nih.gov/pub/cobalt/executables/LATEST/*tar.gz .
rm -rf ftp.ncbi.nlm.nih.gov
tar -xzvf *x64-linux.tar.gz
mv ncbi-cobalt-* cobalt
rm *tar.gz


#Install Megadock
git clone https://github.com/akiyamalab/MEGADOCK.git $INSTALL_DIR/MEGADOCK
sed -i 's/USE_GPU    := 1/USE_GPU    := 0/' $INSTALL_DIR/MEGADOCK/Makefile
sed -i 's/USE_MPI    := 1/USE_MPI    := 0/' $INSTALL_DIR/MEGADOCK/Makefile
sed -i 's#FFTW_INSTALL_PATH ?= /usr/local#FFTW_INSTALL_PATH ?= /usr/local/${FFTW_TARGET}#' $INSTALL_DIR/MEGADOCK/Makefile
cd $INSTALL_DIR/MEGADOCK
make -j$(nproc)

# Clone the PeSTo repository from GitHub
git clone https://github.com/LBM-EPFL/PeSTo.git $INSTALL_DIR/PeSTo
find $INSTALL_DIR/PeSTo/ -name "*.pdb" -type f -delete
rm -rf $INSTALL_DIR/PeSTo/.git
rm -rf $INSTALL_DIR/PeSTo/masif-site_benchmark

## Checkout HADDOCK 3 and custom cns1.3 files
git clone --recursive https://github.com/haddocking/haddock3.git $INSTALL_DIR/haddock3 
cd $INSTALL_DIR/haddock3 
    ## v3.0.0-beta.5
git checkout 1482d85 
cd $INSTALL_DIR
/bin/bash -c "source activate monviso && cd $INSTALL_DIR/haddock3 && pip install -r  requirements.txt && python setup.py install"
git clone --recursive https://github.com/colbyford/HADDOCKer.git /tmp/HADDOCKer
mv /tmp/HADDOCKer/HADDOCK3/cns_solve_1.3_all_intel-mac_linux.tar.gz $INSTALL_DIR

## Note: The custom cns1.3 files are from the HADDOCK3 repo under "varia"
CNS=$INSTALL_DIR/cns_solve 
tar -xzf $INSTALL_DIR/cns_solve_1.3_all_intel-mac_linux.tar.gz -C $INSTALL_DIR/
mv $INSTALL_DIR/cns_solve_1.3/ $CNS 
rm $INSTALL_DIR/cns_solve_1.3_all_intel-mac_linux.tar.gz 
cp $INSTALL_DIR/haddock3/varia/cns1.3/* $CNS/source 
sed -i 's/_CNSsolve_location_/'$INSTALL_DIR'\/cns_solve/g' $CNS/cns_solve_env 
cd $CNS 
make install compiler=gfortran

cd $INSTALL_DIR
#install msms
curl -L 'https://ccsb.scripps.edu/msms/download/933/' --output 'msms_i86_64Linux2_2.6.1.tar.gz'
mkdir -p msms && tar -xvf msms_i86_64Linux2_2.6.1.tar.gz -C msms
rm msms_i86_64Linux2_2.6.1.tar.gz
mv msms/msms.x86_64Linux2.2.6.1 msms/msms

# Download and decompress the UniProt/SwissProt databases into "$INSTALL_DIR"
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz 
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz
gzip -f -d uniprot_sprot.fasta.gz
gzip -f -d uniprot_sprot_varsplic.fasta.gz
#Same for PDB sequences
wget https://files.rcsb.org/pub/pdb/derived_data/pdb_seqres.txt.gz
gzip -f -d pdb_seqres.txt.gz

echo "DB_LOCATION=$INSTALL_DIR" > parameters.txt && \
echo "MEGADOCK_HOME=$INSTALL_DIR/MEGADOCK" >> parameters.txt && \
echo "HDOCKLITE_HOME=$INSTALL_DIR/HDOCKlite-v1.1" >> parameters.txt && \
echo "HADDOCK_HOME=$INSTALL_DIR/haddock3" >> parameters.txt && \
echo "HADDOCK_SELECTION=sasa" >> parameters.txt && \
echo "HADDOCK_CUTOFF=0.9" >> parameters.txt && \
echo "MSMS_HOME=$INSTALL_DIR/msms" >> parameters.txt && \
echo "COBALT_HOME=$INSTALL_DIR/cobalt/bin" >> parameters.txt && \
echo "HMMER_HOME=$INSTALL_DIR/hmmer/bin" >> parameters.txt && \
echo "PESTO_HOME=$INSTALL_DIR/PeSTo" >> parameters.txt && \
echo "MODELLER_EXEC=mod10.5" >> parameters.txt && \
echo "RESOLUTION=4.50" >> parameters.txt && \
echo "SEQID=25" >> parameters.txt && \
echo "HMM_TO_IMPORT=100" >> parameters.txt && \
echo "MODEL_CUTOFF=5" >> parameters.txt && \
echo "PDB_TO_USE=10" >> parameters.txt && \
echo "NUM_OF_MOD_WT=1" >> parameters.txt && \
echo "NUM_OF_MOD_MUT=1" >> parameters.txt && \
echo "W_STRUCT=10" >> parameters.txt && \
echo "W_MUT=10" >> parameters.txt


rm -rf /tmp/HADDOCKer/