#!/bin/bash

java -jar ${HOME}/Applications/womtool-40.jar \
    validate snap-atac-primary.wdl \
    --inputs inputs.json
