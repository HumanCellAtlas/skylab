FROM python:3.6.2
MAINTAINER Jishu Xu <jishuxu@broadinstitute.org>

LABEL software="python 3.6.2"
LABEL description="python 3.6.2 with pysam and scsequtil and more useful modules"

RUN pip3 install scsequtil==0.0.5 \
    scipy \
    numpy \
    matplotlib \
    pandas \
    tables \
    crimson \
    requests

RUN mkdir -p /tools/scripts/
COPY ./parse_picard.py /tools/scripts/
