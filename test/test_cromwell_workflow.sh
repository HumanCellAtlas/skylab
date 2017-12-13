#!/usr/bin/env bash

set -e

CROMWELL_USER=$1
CROMWELL_PASSWORD=$2
CROMWELL_URL=$3
INPUTS_JSON=$4
WDL_FILE=$5
OPTIONS_FILE=$6
DEPENDENCIES_JSON=$7

echo "Starting workflow."
WORKFLOW_HASH=$(cromwell-tools run \
  --username "${CROMWELL_USER}" \
  --password "${CROMWELL_PASSWORD}" \
  --cromwell-url "${CROMWELL_URL}" \
  --inputs-json "${INPUTS_JSON}" \
  --wdl-file "${WDL_FILE}" \
  --options-file "${OPTIONS_FILE}" \
  --dependencies-json "${DEPENDENCIES_JSON}" \
)

echo "Waiting for workflow to complete..."
cromwell-tools wait \
  --username "${CROMWELL_USER}" \
  --password "${CROMWELL_PASSWORD}" \
  --cromwell-url "${CROMWELL_URL}" \
  --workflow-ids "${WORKFLOW_HASH}" \
  --timeout-minutes 60 \
  --poll-interval-seconds 30

echo "Checking workflow completion status."
STATUS_RESULT=$(cromwell-tools status \
  --username "${CROMWELL_USER}" \
  --password "${CROMWELL_PASSWORD}" \
  --cromwell-url "${CROMWELL_URL}" \
  --workflow-ids "${WORKFLOW_HASH}" \
)
echo "${STATUS_RESULT}"

# get the status alone; cromwell tools is a bit more verbose
STATUS=$(echo "${STATUS_RESULT}" | awk '{print $NF}')

if [ "${STATUS}" == "Succeeded" ]; then
  echo "Test Successful. Workflow outputs matched expected values."
  exit 0
else
  echo "Test Failed. Workflow ${WORKFLOW_HASH} did not match expected outputs. See logs for listed workflow for more details."
  exit 1
fi
