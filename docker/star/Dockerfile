FROM ubuntu:17.04
MAINTAINER Jishu Xu <jishuxu@broadinstitute.org>

LABEL software="STAR"
LABEL version="2.5.3a"
LABEL description="RNA-seq aligner"
LABEL website="https://github.com/alexdobin/STAR"
RUN mkdir build 
WORKDIR build 
# install additional python packages
#Install wget, unzip
RUN apt update && apt install -y liblzma-dev libbz2-dev cmake automake curl libboost-all-dev wget build-essential gcc-multilib zlib1g-dev libxml2-dev libncurses5-dev

RUN wget https://github.com/alexdobin/STAR/archive/2.5.3a.tar.gz && \
	tar -xf 2.5.3a.tar.gz
WORKDIR STAR-2.5.3a/bin/Linux_x86_64_static
RUN cp STAR /usr/local/bin

##Install subread 
WORKDIR /build
RUN wget https://downloads.sourceforge.net/project/subread/subread-1.5.3/subread-1.5.3-source.tar.gz
RUN tar -xzvf subread-1.5.3-source.tar.gz
WORKDIR subread-1.5.3-source/src
RUN make -f Makefile.Linux
RUN cp -r /build/subread-1.5.3-source/bin/* /usr/local/bin
RUN ls /build/subread-1.5.3-source/bin/
WORKDIR /


RUN rm -rf /build
