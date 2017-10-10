FROM python:2.7

MAINTAINER Dave Shiga <dshiga@broadinstitute.org> 
LABEL software="Python"
LABEL description="Python with some useful libraries"

RUN mkdir /tools

COPY tools /tools

WORKDIR /tools

RUN cd secondary_analysis && pip install .
