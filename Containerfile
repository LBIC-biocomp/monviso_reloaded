# Use Ubuntu as the base image
FROM ubuntu:latest

# Define a build argument for the HDocklite files
ARG HDOCKLITE_URL
ENV HDOCKLITE_URL=$HDOCKLITE_URL
# Fail the build if the archive URL or path is not provided
RUN if [ -z "$HDOCKLITE_URL" ]; then \
        echo "Error: HDOCKLITE_URL not provided"; \
        exit 1; \
    fi

# Install necessary packages
RUN apt-get update && apt-get install -y wget build-essential git curl make g++ \
    nano cmake gfortran-9 flex csh
RUN apt-get clean 
RUN ln -s /usr/bin/gfortran-9 /usr/bin/gfortran
RUN mkdir /Monviso

# Download and install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh && \
    bash miniconda.sh -b -p /miniconda && \
    rm miniconda.sh

# Add Miniconda to PATH
ENV PATH=/miniconda/bin:${PATH}

# Create a conda environment and install Python 3.9
RUN conda create -n myenv python=3.9 -y
RUN echo "source activate myenv" > ~/.bashrc
ENV PATH /miniconda/envs/myenv/bin:$PATH

# Install Modeller
RUN conda install -c salilab modeller -y

# Define a build argument for MODELLER_LICENSE
ARG MODELLER_LICENSE

# Replace the placeholder in the Modeller configuration file with the license key
RUN sed -i "s/XXXX/${MODELLER_LICENSE}/" /miniconda/lib/modeller-*/modlib/modeller/config.py


# Download and extract Cobalt
RUN wget -r ftp://ftp.ncbi.nlm.nih.gov/pub/cobalt/executables/LATEST/*x64-linux.tar.gz && \
    mv ftp.ncbi.nlm.nih.gov/pub/cobalt/executables/LATEST/*x64-linux.tar.gz . && \
    rm -r ftp.ncbi.nlm.nih.gov && \
    tar -xzvf ncbi-cobalt-*-linux.tar.gz && \
    rm *.tar.gz
RUN mv /ncbi-cobalt-3.0.0 /Monviso/cobalt/

# Install HMMER
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz && \
    tar xvzf hmmer.tar.gz && \
    rm hmmer.tar.gz && \
    hmmer_folder=$(ls | grep hmmer) && \
    cd $hmmer_folder && \
    ./configure --prefix=/Monviso/hmmer/ && \
    make && \
    make install && \
    cd / && \
    rm -r $hmmer_folder

#Install FFTW
ENV FFTW_TARGET fftw-3.3.8

RUN wget -P /tmp http://www.fftw.org/${FFTW_TARGET}.tar.gz && \
    tar xzf /tmp/${FFTW_TARGET}.tar.gz -C /tmp && \
    rm -rf /tmp/${FFTW_TARGET}.tar.gz && \
    cd /tmp/${FFTW_TARGET} && \
    ./configure --enable-float --enable-sse2 --prefix=/usr/local/${FFTW_TARGET} && \
    make -j$(nproc) && \
    make install

#Install Megadock
RUN git clone https://github.com/akiyamalab/MEGADOCK.git /Monviso/MEGADOCK
RUN sed -i 's/USE_GPU    := 1/USE_GPU    := 0/' /Monviso/MEGADOCK/Makefile
RUN sed -i 's/USE_MPI    := 1/USE_MPI    := 0/' /Monviso/MEGADOCK/Makefile
RUN sed -i 's#FFTW_INSTALL_PATH ?= /usr/local#FFTW_INSTALL_PATH ?= /usr/local/${FFTW_TARGET}#' /Monviso/MEGADOCK/Makefile
WORKDIR /Monviso/MEGADOCK
RUN make -j$(nproc)

# Clone the PeSTo repository from GitHub
RUN git clone https://github.com/LBM-EPFL/PeSTo.git /Monviso/PeSTo
RUN find /Monviso/PeSTo/ -name "*.pdb" -type f -delete
RUN rm -r /Monviso/PeSTo/.git
RUN rm -r /Monviso/PeSTo/masif-site_benchmark


## Checkout HADDOCK 3 and custom cns1.3 files
RUN mkdir /Monviso/haddock3

RUN git clone --recursive https://github.com/haddocking/haddock3.git /Monviso/haddock3 && \
    cd /Monviso/haddock3 && \
    ## v3.0.0-beta.5
    git checkout 1482d85 && \
    cd /Monviso
RUN git clone --recursive https://github.com/colbyford/HADDOCKer.git /tmp/HADDOCKer
RUN mv /tmp/HADDOCKer/HADDOCK3/cns_solve_1.3_all_intel-mac_linux.tar.gz /Monviso

## Note: The custom cns1.3 files are from the HADDOCK3 repo under "varia"
RUN export CNS=/Monviso/cns_solve && \
    tar -xzf /Monviso/cns_solve_1.3_all_intel-mac_linux.tar.gz -C /Monviso/&& \
    mv /Monviso/cns_solve_1.3/ $CNS && \
    rm /Monviso/cns_solve_1.3_all_intel-mac_linux.tar.gz && \
    cp /Monviso/haddock3/varia/cns1.3/* $CNS/source && \
    sed -i 's/_CNSsolve_location_/\/Monviso\/cns_solve/g' $CNS/cns_solve_env && \
    chmod -R 777 /Monviso && \
    cd $CNS && \
    make install compiler=gfortran

# Set the working directory
WORKDIR /Monviso

# Download and decompress the UniProt/SwissProt databases into /Monviso
RUN wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && \
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz && \
    gzip -d uniprot_sprot.fasta.gz && \
    gzip -d uniprot_sprot_varsplic.fasta.gz

# Install Monviso using pip in the myenv environment
#RUN git clone https://github.com/LBIC-biocomp/monviso_reloaded
#RUN cd monviso_reloaded
#RUN /bin/bash -c "source activate myenv && pip install -e ."
#RUN cd ..

#Install msms
RUN curl -L 'https://ccsb.scripps.edu/msms/download/933/' --output 'msms_i86_64Linux2_2.6.1.tar.gz'
RUN mkdir -p msms && tar -xvf msms_i86_64Linux2_2.6.1.tar.gz -C msms &&\
rm msms_i86_64Linux2_2.6.1.tar.gz &&\
mv msms/msms.x86_64Linux2.2.6.1 msms/msms

# Copy the HDockLite file and extract it
COPY $HDOCKLITE_URL /Monviso/HDOCKlite.tar.gz
RUN tar -xzf /Monviso/HDOCKlite.tar.gz -C /Monviso; 


RUN echo "DB_LOCATION=/Monviso" > parameters.txt && \
    echo "MEGADOCK_HOME=/Monviso/MEGADOCK" >> parameters.txt && \
    echo "MSMS_HOME=/Monviso/msms" >> parameters.txt && \
    echo "COBALT_HOME=/Monviso/cobalt/bin" >> parameters.txt && \
    echo "HMMER_HOME=/Monviso/hmmer/bin" >> parameters.txt && \
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

RUN rm -rf /tmp/*
RUN rm /Monviso/HDOCKlite.tar.gz