#!/bin/bash

usage()
{
cat << EOF
USAGE: `basename $0` [options]
    -i  inputs.json (e.g. snap-atac.inputs.mouse.json)
EOF
}

while getopts "i:h" OPTION
do
    case $OPTION in
        i) inputs_json=$OPTARG ;;
        h) usage; exit 1 ;;
        *) usage; exit 1 ;;
    esac
done

if [ -z "$inputs_json" ] || [ ! -r "$inputs_json" ]
then
    usage
    exit 1
fi

cromwell-tools submit \
    --url "https://cromwell.caas-prod.broadinstitute.org" \
    --collection-name "pipeline-surge" \
    --service-account-key "${HOME}/pipelines_surge_credentials.json" \
    --wdl snap-atac-primary.wdl \
    --inputs-files ${inputs_json} \
    --options-file snap-atac.options.json
