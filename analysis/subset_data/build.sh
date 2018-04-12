#!/bin/bash

PICARD_JAR=/Users/jsoto/IdeaProjects/picard/build/libs/picard.jar
mkdir docker

cd docker
cp $PICARD_JAR .

echo "FROM java:openjdk-8-jre
MAINTAINER DSDE <dsde@broadinstitute.org>

ENV TERM=xterm-256color

# Change working directory to /usr/gitc/
WORKDIR /root/

# Install python2, python3, and R
RUN apt-get update && \
	apt-get upgrade -y

# Copy everything in working dir outside container to /usr/gitc
COPY . .
" > Dockerfile

# Tagged with the picard release version the jars and binaries were taken from.
docker build -t jsotobroad/gold_data:1.0 .
docker push jsotobroad/gold_data:1.0
