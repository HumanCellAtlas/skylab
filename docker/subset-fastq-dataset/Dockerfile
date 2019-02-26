FROM python:3.6.2

LABEL maintainer="Mint Team <mintteam@broadinstitute.org>" \
      software="python 3.6.2" \
      description="python 3.6.2 with pysam, sctools, requests, and a basic science stack"

RUN apt update && apt install -y \
    samtools

RUN pip3 install \
    Click==7.0 \
    numpy==1.16.1 \
    pysam==0.15.2

RUN mkdir /tools

COPY filterFastqByReadName.py /tools/
COPY prepSubset.sh /tools/

ENV PATH=/tools/:${PATH}

## Note that this string is specific this linux base
ENV CLOUD_SDK_REPO="cloud-sdk-jessie"

RUN echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" |  tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
RUN curl https://packages.cloud.google.com/apt/doc/apt-key.gpg |  apt-key add -
RUN apt update
RUN apt install -y google-cloud-sdk

