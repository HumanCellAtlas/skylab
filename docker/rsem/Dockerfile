FROM centos:centos7

LABEL maintainer="Jishu Xu <jishuxu@broadinstitute.org>" \
  software="RSEM" \
  version="v1.3.0" \
  Description="This image is used for running RSEM RNA-Seq analysis tools " \
  website="https://github.com/deweylab/RSEM"

WORKDIR /usr/local/
RUN yum groupinstall -y 'Development Tools' 
RUN yum install -y perl-Env 
RUN yum install -y zlib-devel
RUN git clone --depth 1 --branch v1.3.0 https://github.com/deweylab/RSEM.git
WORKDIR /usr/local/RSEM
RUN make 
ENV PATH /usr/local/RSEM:$PATH
