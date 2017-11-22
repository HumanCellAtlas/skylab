FROM java:openjdk-8-jre

LABEL maintainer="Jishu Xu <jishuxu@broadinstitute.org>" \
  software="Picard" \
  version="2.10.10" \
  description="A set of command line tools (in Java) for manipulating high-throughput sequencing (HTS) data and formats such as SAM/BAM/CRAM and VCF." \
  website="http://broadinstitute.github.io/picard"

ENV version 2.10.10

# Please follow the below instructions to invoke picard when you are using this docker image:
# java jvm-args -jar /usr/picard/picard.jar PicardToolName OPTION1=value1 OPTION2=value2...

WORKDIR /usr/picard
ADD https://github.com/broadinstitute/picard/releases/download/${version}/picard.jar ./picard.jar
