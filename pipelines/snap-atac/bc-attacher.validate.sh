#!/bin/bash

java -jar ${HOME}/Applications/womtool-40.jar \
    validate bc-attacher.wdl \
    --inputs bc-attacher.inputs.json
