FROM ubuntu:16.04
LABEL maintainer=" Jishu Xu <jishuxu@broadinstitute.org> " \
      software="rsem" \
      version="1.3.0" \
      description="RNA-seq gene expression quantification tools" \
      website="https://deweylab.github.io/RSEM/"

WORKDIR /usr/local/

RUN apt-get update
# Install compiler, perl , R and stuff
RUN apt-get update --fix-missing && \
  apt-get install -y \
  build-essential \
  gcc-multilib \
  apt-utils \
  zlib1g-dev \
  libxml2-dev \
  curl \
  wget \
  git \
  perl \
  perl-base \
  libbz2-dev \
  cmake automake \
  libboost-all-dev \
  libncurses5-dev \
  r-base \
  r-base-core \
  r-base-dev \
  bowtie \
  bowtie2 

# Compile and install STAR
WORKDIR /usr/local
RUN wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz && \
    tar -xf 2.5.3a.tar.gz
WORKDIR STAR-2.5.3a/bin/Linux_x86_64_static
RUN cp STAR /usr/local/bin

# Install RSEM 
WORKDIR /usr/local/
RUN git clone https://github.com/deweylab/RSEM.git
WORKDIR /usr/local/RSEM
RUN git checkout v1.3.0
RUN make 
RUN make ebseq
ENV PATH /usr/local/RSEM:$PATH
