FROM ubuntu:16.04

## Install Variables
ENV pyversion 2.7.13

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

## CellRanger 2.1.1
RUN apt-get -y install git && git clone https://github.com/10XGenomics/cellranger.git
WORKDIR cellranger
RUN make
WORKDIR /root

## Martian
RUN git clone https://github.com/martian-lang/martian.git --recursive
WORKDIR martian
RUN git fetch origin 2.3
RUN git checkout origin/2.3
RUN git submodule update --recursive
RUN make mrc mrf mrg mrp mrs mrt_helper mrstat mrjob
WORKDIR /root

## Set paths
ENV PATH $PATH:/root/cellranger/bin:/root/cellranger/lib/bin:/root/cellranger/tenkit/bin:/root/martian/bin
ENV PYTHONPATH $PYTHONPATH:/root/cellranger/lib/python:/root/cellranger/tenkit/lib/python:/root/martian/adapters/python
ENV MROPATH $MROPATH:/root/cellranger/mro:/root/cellranger/tenkit/mro
ENV RUST_SRC_PATH $RUST_SRC_PATH:/root/cellranger/lib/rust
ENV _TENX_LD_LIBRARY_PATH tenx_path

RUN apt-get update && apt-get -y install python-pip liblzma-dev liblz4-tool libopenblas-dev libatlas-base-dev python-tables cython
COPY requirements.txt /opt/requirements.txt
RUN pip install -r /opt/requirements.txt

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
 && tar xjvf samtools-1.9.tar.bz2 \
 && rm samtools-1.9.tar.bz2 \
 && cd samtools-1.9 \
 && ./configure --prefix=/usr \
 && make \
 && make install \
 && cd .. \
 && rm -rf samtools-1.9

# Install STAR aligner
RUN wget https://github.com/alexdobin/STAR/archive/2.5.1b.tar.gz \
 && tar xf 2.5.1b.tar.gz \
 && rm 2.5.1b.tar.gz \
 && cd STAR-2.5.1b \
 && make \
 && mv bin/Linux_x86_64_static/STAR* /usr/bin \
 && cd .. \
 && rm -rf STAR-2.5.1b

# Install tsne python package. pip installing it doesn't work
RUN git clone https://github.com/danielfrg/tsne.git \
 && cd tsne \
 && make install \
 && cd .. \
 && rm -rf tsne

## Default command
CMD ["cellranger", "-h"]
