#!/usr/bin/env bash

set -e

# NOTE: this script will execute from the repository root when called on jenkins

PIPELINE_FOLDER_NAME=$1
WD=$(pwd)


INPUTS_JSON="/working/test/${PIPELINE_FOLDER_NAME}/pr/test_inputs.json"
WDL_FILE="/working/test/${PIPELINE_FOLDER_NAME}/pr/test_${PIPELINE_FOLDER_NAME}_PR.wdl"
DEPENDENCIES_JSON="/working/test/${PIPELINE_FOLDER_NAME}/pr/dependencies.json"

echo "Setting Cromwell environmental variables"

echo ${BROAD_CROMWELL_KEY} > caas-prod.json
CROMWELL_KEY_FILE="caas-prod.json"

OPTIONS_FILE="https://raw.githubusercontent.com/HumanCellAtlas/skylab/master/test/options.json"
CROMWELL_URL="https://cromwell.caas-prod.broadinstitute.org"
COLLECTION="lira-dev"

echo "Running test"
docker run --rm \
  -v ${WD}:/working \
  -w /working \
  --privileged \
  quay.io/broadinstitute/cromwell-tools:v2.1.0 \
  /working/test/test_cromwell_workflow.sh \
    "${CROMWELL_KEY_FILE}" \
    "${CROMWELL_URL}" \
    "${INPUTS_JSON}" \
    "${WDL_FILE}" \
    "${OPTIONS_FILE}" \
    "${DEPENDENCIES_JSON}" \
    "${COLLECTION}"

