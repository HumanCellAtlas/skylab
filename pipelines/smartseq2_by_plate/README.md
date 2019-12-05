# SmartSeq2 scRNA-Seq Data processing Pipeline for Multiple Cells

This pipeline performs the following tasks:
- Runs the SS2 pipeline for each cell
- Merges the resulting Zarr output files
 
# Inputs
 - Inputs are identical for the SS2 pipeline however, for running 
 multiple samples a URL prefix is provided along with an array of filename prefixes 
 under that path. File names are autogerated by appending _1.fastq.gz and (in the case of 
 paired end) _2.fastq.gz at the end of the filenames
