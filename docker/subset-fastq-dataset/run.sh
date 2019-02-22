#!/bin/bash

./prepSubset.sh \
    -a gs://broad-dsde-mint-dev-cromwell-execution/cromwell-executions/Optimus/4a43c049-39bc-45be-a7a4-fac7e4dc6a41/call-StarAlign \
    -f gs://hca-dcp-sc-pipelines-test-data/referenceData/fastq/chemistry_X10_V2/4k_pbmc/fastqs/ \
    -r chr1
