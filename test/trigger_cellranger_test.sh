#!/usr/bin/env bash

# note, this test must be run from the workflow root for the relative paths to work out.
# WDL file should be the KCO cellranger, located here:
# https://portal.firecloud.org/#methods/single-cell-portal/cell-ranger-2-0-2/3/wdl

set -e

export CROMWELL_SECRETS=~/.ssh/mint_cromwell_config.json

TEST_DIR=./test/optimus/scientific_cellranger

cromwell-tools run \
  --wdl-file "${TEST_DIR}"/cellranger.wdl \
  --inputs-json "${TEST_DIR}"/test_inputs.json
