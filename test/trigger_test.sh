#!/usr/bin/env bash

set -e

# NOTE: this script will execute from the repository root when called on jenkins

VAULT_TOKEN=$1
INPUTS_JSON=$2
WDL_FILE=$3
DEPENDENCIES_JSON=$4
WD=$(pwd)

echo "Rendering templates"
bash ${WD}/test/config/render-ctmpls.sh "dev" "${VAULT_TOKEN}"

echo "Setting Cromwell environmental variables"
CROMWELL_SECRETS="${WD}/test/cromwell-secrets.json"
CROMWELL_USER=$(cat "${CROMWELL_SECRETS}" | jq -r .cromwell_user)
CROMWELL_PASSWORD=$(cat "${CROMWELL_SECRETS}" | jq -r .cromwell_password)

OPTIONS_FILE="https://raw.githubusercontent.com/HumanCellAtlas/skylab/master/test/options.json"
CROMWELL_URL="https://cromwell.mint-dev.broadinstitute.org"

echo "Running test"
docker run --rm \
  -v ${WD}:/working \
  -w /working \
  humancellatlas/cromwell-tools:1.0.1 \
  /working/test/test_cromwell_workflow.sh \
    "${CROMWELL_USER}" \
    "${CROMWELL_PASSWORD}" \
    "${CROMWELL_URL}" \
    "${INPUTS_JSON}" \
    "${WDL_FILE}" \
    "${OPTIONS_FILE}" \
    "${DEPENDENCIES_JSON}"
