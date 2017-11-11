FROM python:2.7

LABEL maintainer="Dave Shiga <dshiga@broadinstitute.org>" \
  software="Python" \
  description="Python2.7 with some useful libraries used for secondary-analysis-python application"

RUN mkdir /tools

COPY tools /tools

WORKDIR /tools

RUN cd secondary_analysis && pip install .

# Install hca-cli for submit.wdl
RUN pip install hca==2.3.0
