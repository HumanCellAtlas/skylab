FROM ubuntu:zesty

LABEL maintainer="Ambrose J. Carr <acarr@broadinstitute.org>" \
  software="kallisto" \
  version="0.43.1" \
  description="kallisto RNA-seq pseudoaligner" \
  website="https://pachterlab.github.io/kallisto"

RUN apt update && apt install -y \
  g++ \
  cmake \
  curl \
  libhdf5-dev \
  samtools

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
