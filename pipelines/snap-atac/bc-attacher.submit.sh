#!/bin/bash

cromwell-tools submit \
    --url "https://cromwell.caas-prod.broadinstitute.org" \
    --collection-name "pipeline-surge" \
    --service-account-key "${HOME}/pipelines_surge_credentials.json" \
    --wdl bc-attacher.wdl \
    --inputs-files bc-attacher.inputs.json \
    --options-file bc-attacher.options.json
