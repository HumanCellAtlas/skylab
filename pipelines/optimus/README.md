# Optimus Pipeline

The Optimus pipeline is a pipeline for processing 3' single-cell expression data generated with the 10X Genomic V2 assay, developed by the Data Coordination Platform (DCP) of the Human Cell Atlas (HCA) Project. 

Optimus is a alignment and transcriptome quantification pipeline. Optimus corrects Cell Barcodes (CBs) and Unique Molecular Identifiers (UMIs), aligns reads to the genome, generates an expression count matrix in a UMI-aware manner, detects empty droplets, calculates summary statistics for genes and cells, and returns outputs in BAM and Zarr file formats. Special care is taken to keep  all reads that may be useful to the downstream user, such as unaligned reads or reads with uncorrectable barcodes. This design provides flexibility to the downstream user and allows for alternative filtering or leveraging the data for novel methodological development.
