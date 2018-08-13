FROM ubuntu:16.04

## Install Variables
ENV pyversion 2.7.13


## Dependencies python 2.7.13
## rust 1.19.0
## clang 4.0.1
## go 1.9

## Python
RUN apt-get update && apt-get -y install build-essential checkinstall libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev tk-dev libc6-dev libbz2-dev wget make && wget https://www.python.org/ftp/python/${pyversion}/Python-${pyversion}.tgz && tar -xvf Python-${pyversion}.tgz
WORKDIR Python-${pyversion}
RUN ./configure
RUN make
WORKDIR /root

## Clang
RUN apt-get -y install clang-4.0 golang-1.9 libz-dev

## RUST start
RUN apt-get -y install curl && curl https://sh.rustup.rs -sSf | sh -s -- -y

## Set Path
ENV PATH $PATH:/Python-${pyversion}:/usr/lib/go-1.9/bin:/root/.cargo/bin

## RUST finish
RUN rustup install 1.19.0 && rustup default 1.19.0

## Node
RUN apt-get update && apt-get install curl
RUN curl -sL https://deb.nodesource.com/setup_10.x | bash -
RUN apt-get install -y build-essential nodejs
RUN npm install npm@latest -g

## CellRanger 2.1.1
RUN apt-get -y install git && git clone https://github.com/10XGenomics/cellranger.git
WORKDIR cellranger
RUN make

## Set path
ENV PATH $PATH:/root/cellranger/bin

## Martian
WORKDIR /root
RUN git clone https://github.com/martian-lang/martian.git --recursive
WORKDIR martian
RUN git fetch origin 2.3
RUN git checkout origin/2.3
RUN git submodule update --recursive
RUN make mrs

## Default command
CMD ["cellranger", "-h"]
