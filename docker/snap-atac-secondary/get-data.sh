#!/bin/bash

# for Step 3
wget -q http://renlab.sdsc.edu/r3fang/share/Fang_2019/MOs_snATAC/genes/promoter.bed -O promoter.bed

# for Step 5
wget -q http://mitra.stanford.edu/kundaje/akundaje/release/blacklists/mm10-mouse/mm10.blacklist.bed.gz -O mm10.blacklist.bed.gz

# for Step 14
wget -q http://renlab.sdsc.edu/r3fang/share/Fang_2019/MOs_snATAC/genes/gencode.vM16.gene.bed -O gencode.vM16.gene.bed
