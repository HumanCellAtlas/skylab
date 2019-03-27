# Optimus Pipeline

The Optimus pipeline, developed by the Data Coordination Platform of the Human Cell Atlas (HCA DCP), processes 3' prime single cell gene expression data from the 10X Genomics v2 (and v3)assay. 

Optimus is a quality control, alignment and transcriptome quantification module. It corrects Cell Barcodes (CBs) and Unique Molecular Identifiers (UMIs), aligns reads to the genome, generates an expression count matrix in a UMI-aware manner, detects empty droplets, calculates summary statistics, and returns outputs in BAM and Zarr file formats. Special care is taken to avoid the removal of reads that are not aligned or that do not contain recognizable barcodes. This design (which differs from many pipelines currently available) allows the use of the entire dataset by those who may want to use alternative filtering or leverage the data for methodological development associated with the data processing.
