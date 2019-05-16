#!/bin/bash

usage()
{
cat << EOF
USAGE: `basename $0` [options]
    -k  service account key (e.g. pipelines_surge_credentials.json)
EOF
}

while getopts "k:h" OPTION
do
    case $OPTION in
        k) service_account_key=$OPTARG ;;
        h) usage; exit 1 ;;
        *) usage; exit 1 ;;
    esac
done

if [ -z "$service_account_key" ] || [ ! -r "$service_account_key" ]
then
    usage
    exit 1
fi

cromwell-tools submit \
    --url "https://cromwell.caas-prod.broadinstitute.org" \
    --collection-name "pipeline-surge" \
    --service-account-key "${service_account_key}" \
    --wdl SnapATACSecondaryAnalysis.wdl \
    --inputs-files SnapATACSecondaryAnalysis.inputs.json
