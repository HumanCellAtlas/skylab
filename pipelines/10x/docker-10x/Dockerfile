FROM ubuntu:zesty

USER root

RUN apt-get update \
 && apt-get install -y alien unzip wget \
 && wget https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-19-1-linux.zip \
 && unzip bcl2fastq2*.zip \
 && alien bcl2fastq2*.rpm \
 && dpkg -i bcl2fastq2*.deb \
 && rm bcl2fastq2*.deb bcl2fastq2*.rpm bcl2fastq2*.zip

ADD cellranger-2.0.0.tar.gz /opt/tenx
COPY run_in_10x_env.bash /opt/tenx/run_in_10x_env.bash

RUN /bin/bash -c "source /opt/tenx/cellranger-2.0.0/sourceme.bash && pip install --upgrade pyparsing"
RUN /bin/bash -c "source /opt/tenx/cellranger-2.0.0/sourceme.bash && pip install martian_cli"

RUN apt-get update \
 && apt-get install -y jq
