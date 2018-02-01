FROM ubuntu:16.04

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="kallisto" \
  version="0.43.1" \
  description="kallisto RNA-seq pseudoaligner" \
  website="https://pachterlab.github.io/kallisto"

RUN apt update && apt install -y \
  g++ \
  cmake \
  curl \
  libhdf5-dev

# Install samtools dependencies
RUN apt install -y \
  libncurses5-dev \
  openssl \
  liblzma-dev \
  libbz2-dev \
  bzip2 \
  libcurl4-openssl-dev \
  libssl-dev

RUN mv /usr/include/hdf5/serial/* /usr/include/

# download, make and install the kallisto binary
RUN curl -L -o kallisto-0.43.1.tar.gz \
  https://github.com/pachterlab/kallisto/archive/v0.43.1.tar.gz && \
  tar -xzf kallisto-0.43.1.tar.gz && \
  cd kallisto-0.43.1/ && \
  mkdir build && \
  cd build && \
  cmake .. && \
  make && \
  make install && \
  cd / && \
  rm -r kallisto-0.43.1 && \
  rm kallisto-0.43.1.tar.gz

# Install samtools
WORKDIR /usr/local/samtools
ADD https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 .

RUN tar -xvf samtools-1.6.tar.bz2 && \
    rm samtools-1.6.tar.bz2 && \
    cd samtools-1.6 && \
    ./configure --prefix=/usr && \
    make && \
    make install
