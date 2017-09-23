FROM python:2.7

MAINTAINER Dave Shiga <dshiga@broadinstitute.org> 
LABEL software="Python"
LABEL description="Python with some useful libraries"

RUN pip install \
    requests \
    boto3

RUN mkdir /tools

COPY tools /tools

RUN pip install -e '/tools/submit'
