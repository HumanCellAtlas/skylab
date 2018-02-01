FROM ubuntu:16.04

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="STAR" \
  version="2.5.3a-40ead6e" \
  description="RNA-seq aligner, patch of version 2.5.3a to allow uBam alignment" \
  website="https://github.com/alexdobin/STAR"

RUN apt update && apt install -y \
  g++ \
  cmake \
  curl \
  git

# Install samtools dependencies
RUN apt install -y \
  libncurses5-dev \
  openssl \
  liblzma-dev \
  libbz2-dev \
  bzip2 \
  libcurl4-openssl-dev \
  libssl-dev

RUN git clone https://github.com/alexdobin/STAR.git && \
  cd STAR && \
  git reset --hard 40ead6efcd9ee2ebf77a510283b69b1a7aad1183 && \
  cp ./bin/Linux_x86_64_static/STAR /usr/local/bin && \
  cp ./LICENSE .. && \
  cd .. && \
  rm -rf STAR

# Install samtools
WORKDIR /usr/local/samtools
ADD https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 .

RUN tar -xvf samtools-1.6.tar.bz2 && \
    rm samtools-1.6.tar.bz2 && \
    cd samtools-1.6 && \
    ./configure --prefix=/usr && \
    make && \
    make install
