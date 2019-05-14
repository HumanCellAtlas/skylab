# Start with an ubuntu system and update it
FROM ubuntu:18.10

# Image label
LABEL maintainer="Nick Barkas <nbarkas@broadinstititute.org>" \
  software="Snaptools" \
  version="0.0.1" \
  description="Image with R supporting snap-atac" \
  website="https://github.com/r3fang/SnapATAC/"

ENV DEBIAN_FRONTEND=noninteractive

# Enable source repositories to install deps for R and update the apt-get list
RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list && apt-get update

RUN apt-get install -y r-base libcurl4-gnutls-dev libxml2-dev libssl-dev git

# install R dependencies and clean up
COPY installRDeps.R /root/
RUN Rscript /root/installRDeps.R

RUN git clone https://github.com/r3fang/SnapATAC.git && cd SnapATAC && R CMD INSTALL .