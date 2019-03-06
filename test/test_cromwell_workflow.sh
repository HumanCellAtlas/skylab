#!/usr/bin/env bash

set -e

CROMWELL_USER=$1
CROMWELL_PASSWORD=$2
CROMWELL_URL=$3
INPUTS_JSON=$4
WDL_FILE=$5
OPTIONS_FILE=$6
DEPENDENCIES_JSON=$7
# Read list of dependency files into an array to pass to cromwell-tools
DEPENDENCIES_LIST=$(cat ${DEPENDENCIES_JSON} | python3 -c "import json,sys;obj=json.load(sys.stdin);print(' '.join(obj.values()));")

echo "Starting workflow."
WORKFLOW_HASH=$(cromwell-tools submit \
  --username "${CROMWELL_USER}" \
  --password "${CROMWELL_PASSWORD}" \
  --url "${CROMWELL_URL}" \
  --inputs_files "${INPUTS_JSON}" \
  --wdl-file "${WDL_FILE}" \
  --options-file "${OPTIONS_FILE}" \
  --deps-file ${DEPENDENCIES_LIST[@]} \
)

# Get workflow id from cromwell-tools response: {"id": "XXXXX", "status": "Submitted"}
WORKFLOW_ID=$(echo ${WORKFLOW_HASH} | python3 -c "import json,sys;obj=json.load(sys.stdin);print(obj['id']);")

echo "Waiting for workflow ${WORKFLOW_ID} to complete..."
cromwell-tools wait "${WORKFLOW_ID}" \
  --username "${CROMWELL_USER}" \
  --password "${CROMWELL_PASSWORD}" \
  --url "${CROMWELL_URL}" \
  --timeout-minutes 60 \
  --poll-interval-seconds 30

echo "Checking workflow completion status."
STATUS_RESULT=$(cromwell-tools status \
  --username "${CROMWELL_USER}" \
  --password "${CROMWELL_PASSWORD}" \
  --url "${CROMWELL_URL}" \
  --uuid "${WORKFLOW_ID}" \
)

# get the status alone; cromwell tools is a bit more verbose
STATUS=$(echo ${STATUS_RESULT} | python3 -c "import json,sys;obj=json.load(sys.stdin);print(obj['status']);")

if [ "${STATUS}" == "Succeeded" ]; then
  echo "Test Successful. Workflow outputs matched expected values."
  exit 0
else
  echo "Test Failed. Workflow ${WORKFLOW_HASH} did not match expected outputs. See logs for listed workflow for more details."
  exit 1
fi
