#!/bin/bash

docker build . -t quay.io/humancellatlas/snaptools:0.0.1
docker push quay.io/humancellatlas/snaptools:0.0.1
