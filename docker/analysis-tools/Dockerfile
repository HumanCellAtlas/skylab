FROM openjdk:8-jre

LABEL maintainer="Jishu Xu <jishuxu@broadinstitute.org>" \
    software="Analysis-tools with Picard-2.10.10, Python-3.5.3 and R-3.3.3" \
    description="A generic toolset for doing large-scale analysis easily with Google Cloud Buckets."

# Install Picard

# Please follow the below instructions to invoke picard when you are using this docker image:
# java jvm-args -jar /usr/picard/picard.jar PicardToolName OPTION1=value1 OPTION2=value2...
ENV picard_version 2.10.10
WORKDIR /usr/picard
ADD https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar ./picard.jar

# Install Python3
ENV PATH /usr/local/bin:$PATH
RUN apt-get update \
  && apt-get install -y python3-pip python3-dev \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 install --upgrade pip

# Install R and other dependencies
RUN apt update && apt install -y \
  	build-essential \
	cmake automake \
  	curl \
  	gcc-multilib \
  	git \
  	libcurl4-openssl-dev \
    libssl-dev \
  	libboost-all-dev \
  	libncurses5-dev \
	libxml2-dev \
  	libncurses5-dev \
  	libboost-all-dev \
	libbz2-dev \
	liblzma-dev \
	lsb-release \
	r-base \
  	r-base-core \
  	r-base-dev \
	samtools \
	sudo \
	wget \
	zlib1g-dev

# Install python packages
RUN pip3 install \
    crimson==0.3.0 \
    HTSeq==0.9.0 \
    matplotlib==2.1.0 \
    numpy==1.12.0 \
    pandas==0.20.3 \
    pysam==0.12.0.1 \
    requests==2.18.4 \
    scipy==0.18.1 \
    sctools==0.1.4 \
    tables==3.4.2 \
    google-cloud-storage \
    git+git://github.com/HumanCellAtlas/pipeline-tools.git

# Fix cannot import name 'opentype' error
RUN pip3 install --upgrade google-auth-oauthlib

# Install R packages
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('reshape')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('googleCloudStorageR')"

# Mount Google Cloud secret key
ADD broad-dsde-mint-dev-analysis-tools.json /usr/secrets/broad-dsde-mint-dev-analysis-tools.json

# Set the environment variable for Python
ENV GOOGLE_APPLICATION_CREDENTIALS=/usr/secrets/broad-dsde-mint-dev-analysis-tools.json

# Set the environment variable for R
ENV GCS_AUTH_FILE=$GOOGLE_APPLICATION_CREDENTIALS

# Install gcloud components
RUN curl -sSL https://sdk.cloud.google.com > /tmp/gcl && bash /tmp/gcl --install-dir=/usr/gcloud --disable-prompts

# Configure gcloud
ENV PATH $PATH:/usr/gcloud/google-cloud-sdk/bin
RUN gcloud auth activate-service-account --key-file=$GCS_AUTH_FILE

# Add Cromwell Credentials
ADD broad-dsde-mint-dev-cromwell.json /usr/secrets/broad-dsde-mint-dev-cromwell.json

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz \
	&& tar -zxvf sratoolkit.2.8.2-1-ubuntu64.tar.gz \
	&& cp -r  sratoolkit.2.8.2-1-ubuntu64/ /usr/local/
ENV PATH /usr/local/sratoolkit.2.8.2-1-ubuntu64/bin:$PATH

RUN pip3 install ipython