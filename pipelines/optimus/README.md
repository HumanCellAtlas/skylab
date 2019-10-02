# Optimus Pipeline Overview
![Diagram](Optimus_diagram.png)

## Introduction to the Optimus Workflow

Optimus is a pipeline developed by the Data Coordination Platform (DCP) of the Human Cell Atlas (HCA) Project that supports processing of any 3' single-cell expression data generated with the [10X Genomic V2 and V3 assay](https://www.10xgenomics.com/solutions/single-cell/). It is an alignment and transcriptome quantification pipeline that corrects Cell Barcodes (CBs) and Unique Molecular Identifiers (UMIs), aligns reads to the genome, generates an expression count matrix in a UMI-aware manner, detects empty droplets, calculates summary statistics for genes and cells, returns read outputs in BAM format, and returns expression counts in numpy matrix, Zarr, and loom file formats. Special care is taken to keep all reads that may be useful to the downstream user, such as unaligned reads or reads with uncorrectable barcodes. This design provides flexibility to the downstream user and allows for alternative filtering or leveraging the data for novel methodological development.

Optimus has been validated for both [human]((../../benchmarking/optimus/optimus_report.rst)and [mouse] (https://docs.google.com/document/d/1_3oO0ZQSrwEoe6D3GgKdSmAQ9qkzH_7wrE7x6_deL10/edit) data.  

More details can be found in the [Optimus Acceptance Report in the repository](../../benchmarking/optimus/optimus_report.rst), or [in the original file](https://docs.google.com/document/d/158ba_xQM9AYyu8VcLWsIvSoEYps6PQhgddTr9H0BFmY/edit).

### Installation  
The Optimus pipeline code can be downloaded by cloning the github repository [Skylab](https://github.com/HumanCellAtlas/skylab/). For the latest release of Optimus please see the realease tags prefixed with "optimus" [here](https://github.com/HumanCellAtlas/skylab/releases)

### Requirements  
Optimus can be deployed using [Cromwell](https://software.broadinstitute.org/wdl/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. Optimus can also be run in [Terra](https://app.terra.bio/#workspaces/help-gatk/HCA_Optimus_Pipeline), a cloud-based analysis platform. In this featured workspace the user will find the Optimus pipeline, configurations, required reference data and other inputs, and example testing data.

### Running Optimus
* [Optimus.wdl](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) in the pipelines/optimus folder,  of the repository, implements the workflow by importing individual tasks in task based WDLs in skylab/library.

#### Inputs
Example input file locations are specified in a json file, e.g., [here](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/example_test_inputs.json). The following table provides information on specific input values.


| Input Name | Description | Allowed Values |
|------------|-------------|----------------|
| chemistry  | chemistry used | "tenX_v2", "tenX_v3" |

#### Outputs
The outputs from the Optimus pipeline can be identified from the outputs of the individual tasks, e.g. [here](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/example_test_outputs.json)

Following are the the types of files produced from the pipeline.

| Output Name | Filename, if applicable | Output Type |Output Format | Notes/Description | Store in Data Store? | Tool |
| ------ | ------ | ------ | ------ | ------ | ------ | ------ |
| pipeline_version | | Version of the processing pipeline run on this data | String | This is passed from the processing WDL to the adapter pipelines to be put into the metadata in the HCA | Yes, in metadata |Lira |
| bam | merged.bam | aligned bam | bam | coordinate sorted | Yes | A few tools; need to address this provenance |
| matrix | sparse_counts.npz | GenexCell expression matrix | Numpy array | | Yes |sctools |
| matrix_row_index | sparse_counts_row_index.npy | Index of cells in expression matrix | Numpy array index | | Yes | sctools |
| matrix_col_index | sparse_counts_col_index.npy | Index of genes in expression matrix | Numpy array index | | Yes | sctools | 
| cell_metrics | merged-cell-metrics.csv.gz | cell metrics | compressed csv | Matrix of metrics by cells | Yes | sctools |
| gene_metrics | merged-gene-metrics.csv.gz | gene metrics | compressed csv | Matrix of metrics by genes | Yes| sctools |
| cell_calls | empty_drops_result.csv | cell calls | csv | | Yes | emptyDrops |
| zarr_output_files | {unique_id}.zarr!.zattrs | | zarr store? sparse matrix? | | Yes | | 
| loom_output_file | output.loom | Loom | Loom | Loom file with expression data and metadata | N/A | N/A |

### Components of Optimus
The source code is available from [Github](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl), an overview of the pipeline can be found on the [HCA Data Portal](https://prod.data.humancellatlas.org/) and the benchmarking that was performed on the pipeline can be found [here](https://docs.google.com/document/d/158ba_xQM9AYyu8VcLWsIvSoEYps6PQhgddTr9H0BFmY/edit#heading=h.calfpviouwbg). Some of the tasks in Optimus use the [sctools](https://github.com/HumanCellAtlas/sctools) library of utilities for large scale distributed single cell data processing, and [Picard](https://broadinstitute.github.io/picard/) tools, a set of command line tools for manipulating high-throughput sequencing data in formats such as SAM/BAM/CRAM and VCF.
