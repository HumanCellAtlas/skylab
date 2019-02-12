FROM python:3.6.2

LABEL maintainer="Mint Team <mintteam@broadinstitute.org>" \
      software="python 3.6.2" \
      description="python 3.6.2 with pysam, sctools, requests, and a basic science stack"

RUN apt update && apt install -y \
    samtools

RUN pip3 install \
    crimson==0.4.0 \
    HTSeq==0.9.0 \
    matplotlib==2.1.0 \
    numpy==1.12.0 \
    pandas==0.20.3 \
    pysam==0.12.0.1 \
    requests==2.18.4 \
    scipy==0.18.1 \
    sctools==0.1.6 \
    tables==3.4.2 \
    numcodecs==0.5.5 \
    zarr==2.2.0

RUN mkdir /tools
WORKDIR /tools

COPY create_zarr_ss2.py .
COPY create_zarr_optimus.py .
