#!/usr/bin/env bash

env=$1
vault_token=$2
working_dir=${3:-$PWD}
docker_image=broadinstitute/dsde-toolbox:latest

echo $(ls ${working_dir}/test/config)

echo "pulling docker image ${docker_image}"
docker pull "${docker_image}"

echo "Creating Cromwell config"
docker run --rm -v ${working_dir}:/working \
    -e VAULT_TOKEN=${vault_token} \
    -e INPUT_PATH=/working/test/config \
    -e OUT_PATH=/working/test \
    "${docker_image}" render-templates.sh ${env}
