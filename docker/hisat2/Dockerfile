FROM ubuntu:17.04
LABEL MAINTAINER="Jishu Xu <jishuxu@broadinstitute.org>"
LABEL software="HISAT2"
LABEL version="2-2.1.0"
LABEL description="RNA-seq aligner"
LABEL website="https://ccb.jhu.edu/software/hisat2/index.shtml"

RUN mkdir -p /opt/tools/
WORKDIR /opt/tools

RUN \
 apt update && \ 
 apt install -y \
	liblzma-dev \
  libbz2-dev \
	cmake \
	automake \
	curl \
	libboost-all-dev \
	wget \
	build-essential \	
	gcc-multilib \
	zlib1g-dev \
	libxml2-dev \
	libncurses5-dev \
	zip unzip \
	samtools \
	git \
	r-base \
	r-base-core \
	r-base-dev
  
RUN \
  wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-source.zip && \
  unzip hisat2-2.1.0-source.zip && \
  cd hisat2-2.1.0 && \
  make && \
  cp hisat2* /usr/local/bin

COPY  ./make_grch38_snp_tran_*.sh hisat2-2.1.0/
ENV PATH="/opt/tools/hisat2-2.1.0:${PATH}"


WORKDIR /opt/tools/gffread
RUN \
	git clone https://github.com/gpertea/gclib && \
	git clone https://github.com/gpertea/gffread && \
	cd gffread && \
	make && \ 
	cp gffread /usr/local/bin/

WORKDIR /opt/tools
