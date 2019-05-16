#!/bin/bash

java -jar ${HOME}/Applications/womtool-40.jar \
    validate SnapATACSecondaryAnalysis.wdl \
    --inputs SnapATACSecondaryAnalysis.inputs.json
