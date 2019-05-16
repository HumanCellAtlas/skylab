#!/bin/bash

# pre-download the necessary dataset for secondary analysis
# then build them into the docker image
./get-data.sh

# build and push
docker build . -t quay.io/humancellatlas/snap-atac-secondary:0.0.1
docker push quay.io/humancellatlas/snap-atac-secondary:0.0.1
