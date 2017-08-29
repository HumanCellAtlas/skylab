FROM centos:centos7
MAINTAINER Jishu Xu <jishuxu@broadinstitute.org>
LABEL Description="This image is used for running RSEM RNA-Seq analysis tools "
WORKDIR /usr/local/
RUN yum groupinstall -y 'Development Tools' 
RUN yum install -y perl-Env 
RUN yum install -y zlib-devel
RUN git clone https://github.com/deweylab/RSEM.git
WORKDIR /usr/local/RSEM
RUN make 
ENV PATH /usr/local/RSEM:$PATH
