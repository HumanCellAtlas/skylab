FROM ubuntu:16.04

LABEL maintainer="Mint Team <mintteam@broadinstitute.org>" \
      software="gtfToGenePred" \
      description="gtfToGenePred from UCSC for building refflat files for Picard"

RUN apt update && apt install -y \
    wget\
    libkrb5-dev


RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred

RUN cp gtfToGenePred /usr/local/bin

RUN chmod 777 /usr/local/bin/gtfToGenePred