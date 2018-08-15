# SmartSeq2 scRNA-Seq Data processing Pipeline + Pipeline Benchmarking workflow

This pipeline performs the following tasks:
- Aligns reads of each sample to the reference genome or transcriptome (using HISAT2)
- Performs quality control on generated BAM files (using Picard)
- Calculates gene expression estimation (using RSEM)
- Group and output QCs, gene/isoform epxression levels
- Benchmark the pipeline in this workflow with existing pipelines. 
 
# Inputs
 - Sequencing data. Inputs are set of paired fastq files. The testing datadat are located in a public bucket.
 - Metadata. A metafile include information about samples, such as cell labels, library types.   
-  References. reference and index are required by scientific pipeline and benchmarking tests files.

All inputs are list in `SmartSeq2Plate.json` file
 
# Workflows

The scientific pipeline and benchmarking tests workflows are coded in WDL language. The required dependencies workflows/tasks are listed in `dependency.json` 

# Benchmarking 

The benchmarking tests performance the comparison between two scientific pipelines, the scientific pipeline in this workflow will be benchmarked against with a current existing pipeline. The outcomes of benchmarking tests will be visualized in a set of html files 
