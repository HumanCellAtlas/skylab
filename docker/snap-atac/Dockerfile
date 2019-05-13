FROM ubuntu:18.04

# Install system requirements
RUN apt update  && \
    apt install -y git python python-pip zlib1g zlib1g-dev wget

# Install SnapTools
RUN mkdir /install && \
    cd /install && \
    git clone https://github.com/r3fang/SnapTools.git && \
    cd SnapTools && \
    git reset --hard 23171b85b175225ca7af0c2c3104159cf6387bbd && \
    pip install .

RUN cd /install && \
    wget -O "bwa-0.7.17.tar.bz2" "https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.17.tar.bz2/download" && \
    tar xvjf bwa-0.7.17.tar.bz2 && \
    cd bwa-0.7.17 && \
    make && \
    mkdir /tools/ && \
    cp bwa /tools/

ENV PATH /tools/:$PATH