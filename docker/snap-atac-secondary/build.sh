#!/bin/bash

docker build . -t quay.io/humancellatlas/snap-atac-secondary:0.0.1
docker push quay.io/humancellatlas/snap-atac-secondary:0.0.1
