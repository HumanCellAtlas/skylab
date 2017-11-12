#!/usr/bin/env bash

# Usage:
# bash deploy_pipelines.sh env infra_dir 10x_dir ss2_dir
#
# Where:
#   -env is the environment, used in bucket name, e.g. dev or staging
#   -infra_dir is the name of the subdirectory where wrapper wdls will be deployed
#   -10x_dir is the name of the subdirectory where the 10x wdls will be deployed 
#   -ss2_dir is the name of the subdirectory where the smartseq2 wdls will be deployed
#
# Example:
# bash deploy_pipelines.sh dev ds_fix_submit_123 1.0.1 ah37fh9
#
# would deploy things as follows:
#
#  gs://broad-dsde-mint-dev-pipelines/10x/1.0.1/analysis.wdl
#  gs://broad-dsde-mint-dev-pipelines/10x/1.0.1/infrastructure/ds_fix_submit_123/submit.wdl
#  gs://broad-dsde-mint-dev-pipelines/10x/1.0.1/infrastructure/ds_fix_submit_123/wrapper.wdl
#  gs://broad-dsde-mint-dev-pipelines/10x/1.0.1/infrastructure/ds_fix_submit_123/dependencies.zip
#  gs://broad-dsde-mint-dev-pipelines/10x/1.0.1/infrastructure/ds_fix_submit_123/static_inputs.json
#  gs://broad-dsde-mint-dev-pipelines/10x/1.0.1/infrastructure/ds_fix_submit_123/options.json
#
#  gs://broad-dsde-mint-dev-pipelines/ss2/ah37fh9/analysis.wdl
#  gs://broad-dsde-mint-dev-pipelines/ss2/ah37fh9/infrastructure/ds_fix_submit_123/submit.wdl
#  gs://broad-dsde-mint-dev-pipelines/ss2/ah37fh9/infrastructure/ds_fix_submit_123/wrapper.wdl
#  gs://broad-dsde-mint-dev-pipelines/ss2/ah37fh9/infrastructure/ds_fix_submit_123/dependencies.zip
#  gs://broad-dsde-mint-dev-pipelines/ss2/ah37fh9/infrastructure/ds_fix_submit_123/static_inputs.json
#  gs://broad-dsde-mint-dev-pipelines/ss2/ah37fh9/infrastructure/ds_fix_submit_123/options.json

set -e

env=$1
infra_dir=$2
tenx_dir=$3
ss2_dir=$4
bucket=broad-dsde-mint-$env-pipelines

echo $env
echo $infra_dir
echo $tenx_dir
echo $ss2_dir
echo $bucket

# Go to root of skylab repo
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $script_dir
cd ../
skylab_root=$(pwd)

function copy_to_bucket()
{
  bucket=$1
  pipeline_name=$2
  pipeline_version=$3
  infra_dir=$4
  stage_dir=$5
  gsutil cp $stage_dir/analysis.wdl gs://$bucket/$pipeline_name/$pipeline_version/
  gsutil cp $stage_dir/submit.wdl gs://$bucket/$pipeline_name/$pipeline_version/infrastructure/$infra_dir/
  gsutil cp $stage_dir/wrapper.wdl gs://$bucket/$pipeline_name/$pipeline_version/infrastructure/$infra_dir/
  gsutil cp $stage_dir/dependencies.zip gs://$bucket/$pipeline_name/$pipeline_version/infrastructure/$infra_dir/
  gsutil cp $stage_dir/static_inputs.json gs://$bucket/$pipeline_name/$pipeline_version/infrastructure/$infra_dir/
  gsutil cp $stage_dir/options.json gs://$bucket/$pipeline_name/$pipeline_version/infrastructure/$infra_dir/options.json
  gsutil ls gs://broad-dsde-mint-$env-pipelines/$pipeline_name/$pipeline_version/
  gsutil ls gs://broad-dsde-mint-$env-pipelines/$pipeline_name/$pipeline_version/infrastructure/$infra_dir
}

function make_zip_file()
{
  stage_dir=$1
  alt_wdl_name=$2
  cd $stage_dir
  if [ -e dependencies.zip ]; then
    rm dependencies.zip
  fi
  zip dependencies.zip analysis.wdl submit.wdl $alt_wdl_name
  cd -
}


### Deploy 10x ###

# Copy files to stage dir
stage_dir=$skylab_root/10x/hca-dcp/stage
local_infra_dir=$skylab_root/10x/hca-dcp
analysis_dir=$skylab_root/10x/count
mkdir -p $stage_dir
cp $analysis_dir/count.wdl $stage_dir/analysis.wdl
cp $analysis_dir/count.wdl $stage_dir/
cp $skylab_root/hca-dcp/submit.wdl $stage_dir
cp $local_infra_dir/wrapper_10x_count.wdl $stage_dir/wrapper.wdl
cp $local_infra_dir/wrapper_10x_count_example_static.json $stage_dir/static_inputs.json
cp $local_infra_dir/10x_options.json $stage_dir/options.json

# Make zip file
make_zip_file $stage_dir count.wdl

# Copy files to bucket
copy_to_bucket $bucket 10x $tenx_dir $infra_dir $stage_dir


### Deploy Smart-seq2 ###

# Copy files to stage dir
stage_dir=$skylab_root/smartseq2_single_sample/hca-dcp/stage
local_infra_dir=$skylab_root/smartseq2_single_sample/hca-dcp
analysis_dir=$skylab_root/smartseq2_single_sample
mkdir -p $stage_dir
cp $analysis_dir/ss2_single_sample.wdl $stage_dir/analysis.wdl
cp $analysis_dir/ss2_single_sample.wdl $stage_dir/
cp $skylab_root/hca-dcp/submit.wdl $stage_dir
cp $local_infra_dir/wrapper_ss2_single_sample.wdl $stage_dir/wrapper.wdl
cp $local_infra_dir/wrapper_ss2_single_sample_example_static.json $stage_dir/static_inputs.json
cp $local_infra_dir/ss2_options.json $stage_dir/options.json

# Make zip file
make_zip_file $stage_dir ss2_single_sample.wdl

# Copy files to bucket
copy_to_bucket $bucket ss2 $ss2_dir $infra_dir $stage_dir

