#!/bin/bash

java -jar ${HOME}/Applications/womtool-40.jar \
    validate snap-atac-primary.wdl \
    --inputs snap-atac.inputs.mouse.json

java -jar ${HOME}/Applications/womtool-40.jar \
    validate snap-atac-primary.wdl \
    --inputs snap-atac.inputs.human.json
