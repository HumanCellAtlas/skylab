#!/usr/bin/env Rscript

## Default repo
local({r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos=r)
})

## Install DropletUtils
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5")
BiocManager::install("GenomicRanges")
BiocManager::install("edgeR")

install.packages(c('RANN','scales','RColorBrewer','foreach','doParallel','igraph','plyr','bigmemory','raster','irlba','Rtsne','doSNOW','devtools'))
