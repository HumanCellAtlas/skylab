#!/bin/bash

tag=$1

if [ -z $tag ]; then
    echo -e "\nYou must provide a tag"
    echo -e "\nUsage: bash build_docker.sh TAG\n"
    exit 1
fi

docker build -t quay.io/humancellatlas/secondary-analysis-ss2-plate-aggregation:$tag .

echo "You can now push with docker push quay.io/humancellatlas/secondary-analysis-ss2-plate-aggregation:$tag"
