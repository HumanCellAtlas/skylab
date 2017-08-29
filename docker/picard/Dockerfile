MAINTAINER Jishu Xu <jishuxu@broadinstitute.org> 
LABEL software="Picard"
LABEL description="This image is for Picard tools"
FROM java:openjdk-8-jre
ENV version 2.10.10
WORKDIR /usr/picard
ADD https://github.com/broadinstitute/picard/releases/download/${version}/picard.jar
