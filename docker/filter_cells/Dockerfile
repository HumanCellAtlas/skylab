# Start with an ubuntu system and update it
FROM ubuntu:16.04

# Image label
LABEL maintainer="Nick Barkas <nbarkas@broadinstititute.org>" \
  software="Filter Cells R Image" \
  version="0.0.1" \
  description="R image for running filter cells"

# Enable source repositories to install deps for R and update the apt-get list
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list && apt-get update

# Install R dependencies, wget, Python and numpy
RUN apt-get install -y \
    build-essential \
    wget \
    python \
    python-pip \
    && apt-get build-dep -y r-base

RUN pip install numpy

# Download and unzip R Build R and add to PATH
WORKDIR /root/
RUN wget https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz \
    && tar xzf R-3.5.1.tar.gz \
    && cd "R-3.5.1" \
    && ./configure --prefix=/opt/R/3.5.1/ --enable-R-shlib --with-blas --with-lapack \
    && make \
    && make install \
    && cd /root/ \
    && rm -r R-3.5.1 \
    && rm R-3.5.1.tar.gz

ENV PATH="/opt/R/3.5.1/bin/:$PATH"

# Install R dependencies and clean up
COPY installRDeps.R /root/
RUN Rscript /root/installRDeps.R \
    && rm installRDeps.R

# Add EmptyDrop scripts
RUN mkdir /root/tools
COPY filter_cells.R ./tools/

ENV PATH="/root/tools/:$PATH"
