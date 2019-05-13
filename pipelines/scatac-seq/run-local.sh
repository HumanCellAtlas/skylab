#!/bin/bash

java -jar cromwell-40.jar \
    run alignment.wdl --inputs inputs.json
