#!/usr/bin/env bash

tag=$1
if [ -z "$tag" ]; then
  echo "Must pass in tag to use for docker image, e.g. test.sh test123"
  exit 1
fi
cd ../../../
docker build -t humancellatlas/secondary-analysis-python:$tag .
cd -

docker run humancellatlas/secondary-analysis-python:$tag bash -c "python -m unittest discover -v tests"
