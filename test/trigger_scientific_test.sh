#!/usr/bin/env bash

# Note, when using default parameters, this test must be run from the workflow root for the
# relative paths to work out.

set -e
# This secrets file contains the following information on a cromwell instance:
# username
# password
# url (e.g. https://cromwell.my-url.broadinstitute.org)
CROMWELL_SECRETS=${2:-${HOME}/.ssh/mint_cromwell_config.json}
TEST_DIR=${3:-./test/optimus/scientific}

export CROMWELL_SECRETS

cromwell-tools run \
  --wdl-file "${TEST_DIR}"/TestOptimusScientific.wdl \
  --dependencies-json "${TEST_DIR}"/dependencies.json \
  --inputs-json "${TEST_DIR}"/test_inputs.json
