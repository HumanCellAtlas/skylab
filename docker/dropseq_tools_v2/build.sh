#!/bin/bash

tag=$1

if [ -z $tag ]; then
    echo -e "\nYou must provide a tag"
    echo -e "\nUsage: bash build_docker.sh TAG\n"
    exit 1
fi

image="quay.io/humancellatlas/secondary-analysis-dropseqtools_v2"

docker build -t $image:$tag .

echo "You can now push with docker push $image:$tag"
