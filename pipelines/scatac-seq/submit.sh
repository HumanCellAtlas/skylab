#!/bin/bash

cromwell-tools submit \
    --url "https://cromwell.caas-prod.broadinstitute.org" \
    --collection-name "pipeline-surge" \
    --service-account-key "pipelines_surge_credentials.json" \
    --wdl "snap-atac-primary.wdl" \
    -i "inputs.json" \
    --options-file options.json
