# Use Ubuntu as the base image
FROM ubuntu:latest

# Install necessary packages
RUN apt-get update && apt-get install -y wget build-essential
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

# Set the working directory
WORKDIR /Monviso

# Download and decompress the UniProt/SwissProt databases into /Monviso
RUN wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && \
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz && \
    gzip -d uniprot_sprot.fasta.gz && \
    gzip -d uniprot_sprot_varsplic.fasta.gz

# Install Monviso using pip in the myenv environment
RUN /bin/bash -c "source activate myenv && pip install monviso==0.1.3"
