#!/bin/bash

cromwell-tools submit \
    --url "https://cromwell.caas-prod.broadinstitute.org" \
    --collection-name "pipeline-surge" \
    --service-account-key "${HOME}/pipelines_surge_credentials.json" \
    --wdl snap-atac-primary.wdl \
    --inputs-files snap-atac.inputs.mouse.json \
    --options-file snap-atac.options.json
