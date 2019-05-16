#!/usr/bin/Rscript

## Default repo
local({r <- getOption("repos")
    r["CRAN"] <- "https://cloud.r-project.org"
    options(repos=r)
})

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rhdf5")
BiocManager::install("GenomicRanges")
BiocManager::install("edgeR")

pkgs = c('RANN', 'scales', 'RColorBrewer', 'foreach', 'doParallel', 'igraph', 'plyr', 'bigmemory', 'raster', 'irlba', 'Rtsne', 'doSNOW', 'devtools', 'leiden', 'umap', 'optparse')
install.packages(pkgs, repos='http://cran.us.r-project.org')
