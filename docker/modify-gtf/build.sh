#!/bin/bash

tag=$1

if [ -z $tag ]; then
    echo -e "\nYou must provide a tag"
    echo -e "\nUsage: bash build_docker.sh TAG\n"
    exit 1
fi

docker build -t quay.io/humancellatlas/modify-gtf:$tag .
docker build -t us.gcr.io/broad-gotc-dev/shpilker_modifygtf:$tag .

echo You can now push with
echo docker push quay.io/humancellatlas/modify-gtf:$tag