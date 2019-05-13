#!/bin/bash

#fixme: skip if already running
java -jar cromwell-40.jar server &

cromwell-tools submit \
    --wdl-file alignment.wdl \
    --inputs-files inputs.json \
    --url http://locahost:8000
