#!/usr/bin/env bash

# NOTE: this script will execute from the repository root when called on jenkins

VAULT_TOKEN=$1
INPUTS_JSON=$2
WDL_FILE=$3
DEPENDENCIES_JSON=$4
WD=$(pwd)
echo $(ls)
echo ${WD}

echo "Rendering templates"
bash ${WD}/test/config/render-ctmpls.sh "dev" "${VAULT_TOKEN}"

echo "Setting Cromwell environmental variables"
CROMWELL_SECRETS="${WD}/test/cromwell-secrets.json"
CROMWELL_USER=$(cat "${CROMWELL_SECRETS}" | jq -r .cromwell_user)
CROMWELL_PASSWORD=$(cat "${CROMWELL_SECRETS}" | jq -r .cromwell_password)

OPTIONS_FILE="https://raw.githubusercontent.com/HumanCellAtlas/skylab/ajc-pr-test/test/options.json"
CROMWELL_URL="https://cromwell.mint-dev.broadinstitute.org"

echo "Running test"
docker run --rm \
  -v ${WD}/test:/test \
  humancellatlas/cromwell-tools:1.0.1 \
  /test/test_cromwell_workflow.sh \
    "${CROMWELL_USER}" \
    "${CROMWELL_PASSWORD}" \
    "${CROMWELL_URL}" \
    "${INPUTS_JSON}" \
    "${WDL_FILE}" \
    "${OPTIONS_FILE}" \
    "${DEPENDENCIES_JSON}"
