#!/bin/bash

docker build . -f ./bc-attacher.Dockerfile -t bc-attacher
# docker build . -f ./bc-attacher.Dockerfile -t quay.io/humancellatlas/bc-attacher:0.0.1
# docker push quay.io/humancellatlas/bc-attacher:0.0.1
