#!/bin/bash

usage()
{
cat << EOF
USAGE: `basename $0` [options]
    -i  inputs.json (e.g. snap-atac.inputs.mouse.json)
    -k  service account key (e.g. pipelines_surge_credentials.json)
EOF
}

while getopts "i:k:h" OPTION
do
    case $OPTION in
        i) inputs_json=$OPTARG ;;
        k) service_account_key=$OPTARG ;;
        h) usage; exit 1 ;;
        *) usage; exit 1 ;;
    esac
done

if [ -z "$inputs_json" ] || [ ! -r "$inputs_json" ]
then
    usage
    exit 1
fi

if [ -z "$service_account_key" ] || [ ! -r "$service_account_key" ]
then
    usage
    exit 1
fi

cromwell-tools submit \
    --url "https://cromwell.caas-prod.broadinstitute.org" \
    --collection-name "pipeline-surge" \
    --service-account-key "${service_account_key}" \
    --wdl snap-atac-primary.wdl \
    --inputs-files ${inputs_json} \
    --options-file snap-atac.options.json
