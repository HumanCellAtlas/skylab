#!/bin/bash

java -jar ${HOME}/Applications/womtool-40.jar \
    validate SnapAtacFromCellRanger.wdl \
    --inputs SnapAtacFromCellRanger.inputs.json
