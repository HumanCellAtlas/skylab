#!/usr/bin/env bash

# note, this test must be run from the workflow root for the relative paths to work out.
# WDL file should be the KCO cellranger, located here:
# https://portal.firecloud.org/#methods/single-cell-portal/cell-ranger-2-0-2/3/wdl

set -e

# This secrets file contains the following information on a cromwell instance:
# username
# password
# url (e.g. https://cromwell.my-url.broadinstitute.org)
CROMWELL_SECRETS=${2:-${HOME}/.ssh/mint_cromwell_config.json}
TEST_DIR=${3:-./test/optimus/scientific_cellranger}

export CROMWELL_SECRETS

cromwell-tools run \
  --wdl-file "${TEST_DIR}"/cellranger.wdl \
  --inputs-json "${TEST_DIR}"/test_inputs.json

