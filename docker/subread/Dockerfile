FROM ubuntu:17.04
LABEL MAINTAINER="Jishu Xu <jishuxu@broadinstitute.org>"
LABEL software="subread package"
LABEL version="1.6.0"
LABEL description="RNA-seq high-performance read alignment, quantification and mutation discovery"
LABEL website="http://subread.sourceforge.net/"

# Install compiler 
RUN apt-get update --fix-missing \
	&& apt-get install -y \
	build-essential \
	gcc-multilib \
	apt-utils \
	zlib1g-dev \
	libxml2-dev \
	curl \
	wget \
	libbz2-dev \
	cmake automake \
	libboost-all-dev \
	libncurses5-dev \
	r-base \
	r-base-core \
	r-base-dev \
	libncurses5-dev
        
      
# Install subread 
WORKDIR /usr/local/ 
ENV VERSION="1.6.0"
RUN wget "https://downloads.sourceforge.net/project/subread/subread-${VERSION}/subread-${VERSION}-source.tar.gz"
RUN tar -xzvf subread-${VERSION}-source.tar.gz
WORKDIR /usr/local/subread-${VERSION}-source/src 
RUN make -f Makefile.Linux 
ENV PATH /usr/local/subread-${VERSION}-source/bin/:$PATH
# Cleanup
RUN apt-get clean

