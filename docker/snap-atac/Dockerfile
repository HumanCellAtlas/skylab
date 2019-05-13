FROM ubuntu:18.04

# Install system requirements
RUN apt update  && \
    apt install -y git python python-pip zlib1g zlib1g-dev

# Install SnapTools
RUN mkdir /install && \
    cd /install && \
    git clone https://github.com/r3fang/SnapTools.git && \
    cd SnapTools && \
    git reset --hard 23171b85b175225ca7af0c2c3104159cf6387bbd && \
    pip install .

# This is a separate step becase it will be removed in a later iteration of the pipeline
RUN apt install -y bwa

