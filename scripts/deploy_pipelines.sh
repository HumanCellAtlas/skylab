#!/usr/bin/env bash

set -e

env=$1
infra_version=$2
tenx_version=$3
ss2_version=$4
bucket=broad-dsde-mint-$env-pipelines

echo $env
echo $infra_version
echo $tenx_version
echo $ss2_version
echo $bucket

# Deploy 10x
cd ../10x/hca-dcp
mkdir -p deploy
cd deploy
cp ../../count/count.wdl .
cp ../../../hca-dcp/submit.wdl .
if [ -e dependencies.zip ]; then
  rm dependencies.zip
fi
zip dependencies.zip count.wdl submit.wdl
gsutil cp count.wdl gs://$bucket/10x/$tenx_version/
gsutil cp submit.wdl gs://$bucket/10x/$tenx_version/infrastructure/$infra_version/submit.wdl
gsutil cp ../wrapper_10x_count.wdl gs://$bucket/10x/$tenx_version/infrastructure/$infra_version/wrapper.wdl
gsutil cp ../wrapper_10x_count.wdl gs://$bucket/10x/$tenx_version/infrastructure/$infra_version/wrapper.wdl
gsutil cp dependencies.zip gs://$bucket/10x/$tenx_version/infrastructure/$infra_version/dependencies.zip
gsutil cp ../wrapper_10x_count_example_static.json gs://$bucket/10x/$tenx_version/infrastructure/$infra_version/static_inputs.json
gsutil cp ../10x_options.json gs://$bucket/10x/$tenx_version/infrastructure/$infra_version/options.json
gsutil ls gs://broad-dsde-mint-$env-pipelines/10x/$tenx_version/
gsutil ls gs://broad-dsde-mint-$env-pipelines/10x/$tenx_version/infrastructure/$infra_version
cd ../../../

# TODO: Deploy Smart-seq2
