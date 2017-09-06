FROM ubuntu:17.04
MAINTAINER Ambrose Carr <mail@ambrosejcarr.com>

LABEL software="python3.6"
LABEL description="python3.6 with pysam and scsequtil"

RUN apt-get update && apt-get install -y \
  curl \
  python3.5 \
  python3-pip \
  python3-pysam \
  python3-cairocffi \
  wget

# bind python to python3.5 as well so either can be used
RUN ln -s /usr/bin/python3.5 /usr/bin/python

RUN python3.5 -m pip install scsequtil

WORKDIR /
