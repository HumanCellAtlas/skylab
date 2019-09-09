#MOUSE
cromwell-tools submit --url "https://cromwell.caas-prod.broadinstitute.org" --collection-name "pipeline-surge" --service-account-key /Users/kkonwar/.ssh/pipelines_surge_credentials.json --wdl BuildIndices.wdl -i mouse_inputs.json -l mouse-build-label.json 
