FROM openjdk:8-jre
  
LABEL maintainer="Jishu Xu <jishuxu@broadinstitute.org>" \
  software="DropSeq-tools-HCA" \
  version="1.3.1" \
  description="Dropseq tool is a java program to process Dropseq data " \
  website="https://github.com/broadinstitute/drop-seq-tools-hca"

ENV version 1.3.1
WORKDIR /usr/dropseq-tools-hca/jar
RUN mkdir jar
COPY jar/ jar/

