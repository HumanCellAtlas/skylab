# Start with an ubuntu system and update it
FROM ubuntu:16.04

# Image label
LABEL maintainer="Nick Barkas <nbarkas@broadinstititute.org>" \
  software="DropletUtils" \
  version="1.2.1-0.1.0" \
  description="Bioconductor DropletUtils Package with the suitable version of R" \
  website="https://bioconductor.org/packages/release/bioc/html/DropletUtils.html"

# Enable source repositories to install deps for R
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list
RUN apt-get update

# Build R version 3.5.1 from source
RUN apt-get install -y build-essential

# Install r dependencies
RUN apt-get build-dep -y r-base

# Install wget
RUN apt-get install -y wget

# Download and unzip R
WORKDIR /root/
RUN wget https://cloud.r-project.org/src/base/R-3/R-3.5.1.tar.gz && tar xzf R-3.5.1.tar.gz

# Build R and add to PATH
WORKDIR /root/R-3.5.1
RUN ./configure --prefix=/opt/R/3.5.1/ --enable-R-shlib --with-blas --with-lapack
RUN make && make install
ENV PATH="/opt/R/3.5.1/bin/:$PATH"

# Install DropUtils
WORKDIR /root/
COPY installDropletUtils.R /root/
RUN Rscript /root/installDropletUtils.R

# Install emptyDropsWrapper and npz2rds
RUN apt-get install -y git
RUN git clone https://github.com/HumanCellAtlas/skylab-analysis.git && cd skylab-analysis && git checkout b772460458388c8edaaee596db08829bd1b64cec
ENV PATH="/root/skylab-analysis/tools/emptyDropsWrapper/:/root/skylab-analysis/tools/npz2rds/:$PATH"