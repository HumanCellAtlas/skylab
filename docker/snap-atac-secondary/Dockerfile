FROM ubuntu:18.10

# Image label
LABEL maintainer="Nick Barkas <nbarkas@broadinstititute.org>" \
  contributor="Jaeyoung Chun (chunj@mskcc.org)" \
  software="Snaptools" \
  version="0.0.1" \
  description="SnapATAC secondary analysis" \
  website="https://github.com/r3fang/SnapATAC/"

# non-interactive
ENV DEBIAN_FRONTEND=noninteractive

# enable source repositories to install deps for R and update the apt-get list
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list && apt-get update

RUN apt-get install -y build-essential
RUN apt-get install -y r-base

RUN apt-get install -y libcurl4-openssl-dev libssl-dev
RUN apt-get install -y git

#--> install SnapATAC
COPY install-packages.R /opt/install-packages.R

RUN Rscript /opt/install-packages.R

# fixme: cd /opt && \
RUN git clone https://github.com/r3fang/SnapATAC.git && \
    cd SnapATAC && \
    R CMD INSTALL .
#<--

# fixme: move up later
RUN apt-get install -y python python-pip

#--> install SnapTools
RUN cd /opt && \
    git clone https://github.com/r3fang/SnapTools.git && \
    cd SnapTools && \
    pip install .

# this is a temporary hack
RUN mkdir -p /home/r3fang/anaconda2/bin/ && \
    ln -s `which python` /home/r3fang/anaconda2/bin/python
#<--

# install macs2
RUN pip install --trusted-host pypi.python.org MACS2==2.1.2.1

COPY data /opt/data/
COPY secondary-analysis.R /opt/

# ENTRYPOINT ["Rscript", "/opt/secondary-analysis.R"]
