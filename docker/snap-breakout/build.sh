#!/bin/bash

docker build -t quay.io/humancellatlas/snap-breakout:testing .
docker push quay.io/humancellatlas/snap-breakout:testing
