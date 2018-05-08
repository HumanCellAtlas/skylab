#!/usr/bin/env bash

# note, this test must be run from the workflow root for the relative paths to work out.

set -e

export CROMWELL_SECRETS=~/.ssh/mint_cromwell_config.json

TEST_DIR=./test/optimus/scientific_cellranger

cromwell-tools run \
  --wdl-file "${TEST_DIR}"/cellranger.wdl \
  --inputs-json "${TEST_DIR}"/test_inputs.json
