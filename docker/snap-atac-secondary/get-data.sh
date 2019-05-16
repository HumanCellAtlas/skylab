#!/bin/bash

# pre-download the necessary dataset for secondary analysis
# then build them into the docker image

path_data="./data"
mkdir -p ${path_data}

# for Step 3
wget -q http://renlab.sdsc.edu/r3fang/share/Fang_2019/MOs_snATAC/genes/promoter.bed -O ${path_data}/promoter.bed

# for Step 5
wget -q http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz -O ${path_data}/mm10.blacklist.bed.gz

# for Step 14
wget -q http://renlab.sdsc.edu/r3fang/share/Fang_2019/MOs_snATAC/genes/gencode.vM16.gene.bed -O ${path_data}/gencode.vM16.gene.bed
