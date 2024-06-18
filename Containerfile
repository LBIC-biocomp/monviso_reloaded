# Use Ubuntu as the base image
FROM ubuntu:latest


# Install necessary packages
RUN apt-get update && apt-get install -y wget build-essential git curl
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

# Clone the PeSTo repository from GitHub
RUN git clone https://github.com/LBM-EPFL/PeSTo.git /Monviso/PeSTo
RUN find /Monviso/PeSTo/ -name "*.pdb" -type f -delete
RUN rm -r /Monviso/PeSTo/.git
RUN rm -r /Monviso/PeSTo/masif-site_benchmark

# Set the working directory
WORKDIR /Monviso

# Download and decompress the UniProt/SwissProt databases into /Monviso
RUN wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz && \
    wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot_varsplic.fasta.gz && \
    gzip -d uniprot_sprot.fasta.gz && \
    gzip -d uniprot_sprot_varsplic.fasta.gz

# Install Monviso using pip in the myenv environment
RUN /bin/bash -c "source activate myenv && pip install --default-timeout=200 --retries 5 monviso==0.1.4"

#Install msms
RUN curl --header 'Host: ccsb.scripps.edu'  --header 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' --header 'Accept-Language: en-US,en;q=0.5' --referer 'https://ccsb.scripps.edu/msms/downloads/' --header 'Upgrade-Insecure-Requests: 1' --header 'Sec-Fetch-Dest: document' --header 'Sec-Fetch-Mode: navigate' --header 'Sec-Fetch-Site: same-origin' --header 'Sec-Fetch-User: ?1' 'https://ccsb.scripps.edu/msms/download/933/' --output 'msms_i86_64Linux2_2.6.1.tar.gz' &&\
mkdir -p msms && tar -xvf msms_i86_64Linux2_2.6.1.tar.gz -C msms \
rm msms_i86_64Linux2_2.6.1.tar.gz \
mv msms/msms.x86_64Linux2.2.6.1 msms/msms

RUN echo "DB_LOCATION=/Monviso" > parameters.txt && \
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