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
  samtools \
  sudo \
  wget \
  zlib1g-dev \
  gfortran

## Now install R and littler, and create a link for littler in /usr/local/bin
## Also set a default CRAN repo, and make sure littler knows about it too
RUN echo "deb http://http.debian.net/debian sid main" > /etc/apt/sources.list.d/debian-unstable.list \
  && echo 'APT::Default-Release "testing";' > /etc/apt/apt.conf.d/default
ENV R_BASE_VERSION 3.4.4
RUN apt-get update \
  && apt-get install -t unstable -y --no-install-recommends \
    littler \
    r-cran-littler \
    r-base=${R_BASE_VERSION}-* \
    r-base-dev=${R_BASE_VERSION}-* \
    r-recommended=${R_BASE_VERSION}-* \
  && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
  && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r \
  && ln -s /usr/share/doc/littler/examples/install.r /usr/local/bin/install.r \
  && ln -s /usr/share/doc/littler/examples/install2.r /usr/local/bin/install2.r \
  && ln -s /usr/share/doc/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
  && ln -s /usr/share/doc/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
  && install.r docopt \
  && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
  && rm -rf /var/lib/apt/lists/*

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
# Install gcloud components
RUN curl -sSL https://sdk.cloud.google.com > /tmp/gcl && bash /tmp/gcl --install-dir=/usr/gcloud --disable-prompts
# Configure gcloud
ENV PATH $PATH:/usr/gcloud/google-cloud-sdk/bin


RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2-1/sratoolkit.2.8.2-1-ubuntu64.tar.gz \
  && tar -zxvf sratoolkit.2.8.2-1-ubuntu64.tar.gz \
  && cp -r  sratoolkit.2.8.2-1-ubuntu64/ /usr/local/
ENV PATH /usr/local/sratoolkit.2.8.2-1-ubuntu64/bin:$PATH

RUN pip3 install ipython

# Install R packages
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('reshape')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('googleCloudStorageR')"
RUN Rscript -e "install.packages('gridExtra')"
RUN Rscript -e "install.packages('ggpubr')"
RUN Rscript -e "install.packages('ggpmisc')"
RUN Rscript -e "install.packages('cowplot')"
RUN Rscript -e "install.packages('corrplot')"
RUN Rscript -e "install.packages('ggrepel')"
RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("rtracklayer")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("scran")'
RUN Rscript -e "install.packages('igraph')"
RUN Rscript -e "install.packages('rsvd')"
RUN Rscript -e "install.packages('factoextra')"
RUN Rscript -e "install.packages('fpc')"
RUN Rscript -e "install.packages('NbClust')"
RUN Rscript -e "install.packages('knitr')"
RUN Rscript -e "install.packages('rmarkdown')"
RUN Rscript -e "install.packages('Rtsne')"

# Add benchmarking scripts
RUN apt-get update \
  && apt-get install -t unstable -y --no-install-recommends \
    cabal-install
RUN cabal update && cabal install pandoc
RUN ln -s /root/.cabal/bin/pandoc /usr/local/bin/pandoc
WORKDIR /usr/local/scripts/
ENV PATH /usr/local/scripts/:$PATH
COPY ./*.R /usr/local/scripts/
COPY ./*.py /usr/local/scripts/
