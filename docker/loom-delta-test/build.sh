#!/bin/bash

tag=$1

if [ -z $tag ]; then
    echo -e "\nYou must provide a tag"
    echo -e "\nUsage: bash build_docker.sh TAG\n"
    exit 1
fi

dockerName=loom-delta-test

docker build -t quay.io/humancellatlas/${dockerName}:${tag} .

echo "You can now push with docker push quay.io/humancellatlas/${dockerName}:${tag}"
