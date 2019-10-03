# Optimus Pipeline Overview
![Diagram](Optimus_diagram.png)

## Introduction to the Optimus Workflow

Optimus is a pipeline developed by the Data Coordination Platform (DCP) of the [Human Cell Atlas (HCA) Project](https://data.humancellatlas.org/) that supports processing of any 3' single-cell expression data generated with the [10X Genomic V2 and V3 assay](https://www.10xgenomics.com/solutions/single-cell/). It is an alignment and transcriptome quantification pipeline that corrects Cell Barcodes (CBs), aligns reads to the genome, corrects Unique Molecular Identifiers (UMIs),generates an expression count matrix in a UMI-aware manner, detects empty droplets, calculates summary statistics for genes and cells, returns read outputs in BAM format, and returns expression counts in numpy matrix, Zarr, and loom file formats. Special care is taken to keep all reads that may be useful to the downstream user, such as unaligned reads or reads with uncorrectable barcodes. This design provides flexibility to the downstream user and allows for alternative filtering or leveraging the data for novel methodological development.

Optimus has been validated for analyzing both [human](../../benchmarking/optimus/optimus_report.rst) and [mouse](https://docs.google.com/document/d/1_3oO0ZQSrwEoe6D3GgKdSmAQ9qkzH_7wrE7x6_deL10/edit) data sets. More details about the human validation can be found in the [in the original file](https://docs.google.com/document/d/158ba_xQM9AYyu8VcLWsIvSoEYps6PQhgddTr9H0BFmY/edit).

## Quick Start Table

| Pipeline Features | Description | Source |
|-------------------|---------------------------------------------------------------|-----------------------|
|Assay Type | 10x Single Cell Expression (v2) |[10x Genomics](https://www.10xgenomics.com)
| Overall Workflow  |Quality control module and transcriptome quantification module | Code available from [Github](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) |
| Workflow Language |WDL          |[openWDL](https://github.com/openwdl/wdl)|
| Genomic Reference Sequence|GRCh38 human genome primary sequence and M23 (GRCm38.p6) mouse genome primary sequence|GENCODE [human](https://www.gencodegenes.org/human/release_27.html) and [mouse](https://www.gencodegenes.org/mouse/release_M23.html)|
|Transcriptomic Reference Annotation |V27 GenCode human transcriptome |[GENCODE](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz)|
| Aligner           |STAR       |[Dobin, et al.,2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/)|
| Transcript Quantification |Utilities for processing large-scale single cell datasets |[Sctools](https://github.com/HumanCellAtlas/sctools) |                      
|Data Input File Format |File format in which sequencing data is provided |[FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |                       
|Data Output File Format |File formats in which Optimus output is provided |[BAM](http://samtools.github.io/hts-specs/), [Zarr version 2](https://zarr.readthedocs.io/en/stable/spec/v2.html), Python numpy arrays (internal), [loom](http://loompy.org/) |

# Set-up
## Installation  
The Optimus pipeline code can be downloaded by cloning the github repository [Skylab](https://github.com/HumanCellAtlas/skylab/). For the latest release of Optimus, please see the realease tags prefixed with "optimus" [here](https://github.com/HumanCellAtlas/skylab/releases)

## Requirements  
Optimus can be deployed using [Cromwell](https://software.broadinstitute.org/wdl/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. Optimus can also be run in [Terra](https://app.terra.bio/#workspaces/help-gatk/HCA_Optimus_Pipeline), a cloud-based analysis platform. In this featured workspace the user will find the Optimus pipeline, configurations, required reference data and other inputs, and example testing data.

## Input Files

### Input Data Preparation 

Each 10X v2 3’ sequencing experiment generates triplets of Fastq files:

1. forward reads (R1), containing the unique molecular identifier (UMI) and cell barcode sequences
2. reverse reads (R2), which contain the alignable genomic information from the mRNA transcript 
3. an index fastq file that contains the sample barcodes, when provided by the sequencing facility

Example input file locations are specified in a json file, e.g., [here](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/example_test_inputs.json). The following table provides information on specific input values.


| Input Name | Description | Allowed Values |
|------------|-------------|----------------|
| chemistry  | chemistry used | "tenX_v2", "tenX_v3" |

# Running Optimus

* [Optimus.wdl](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) in the pipelines/optimus folder,  of the repository, implements the workflow by importing individual tasks in task based WDLs in skylab/library.


## Optimus Modules Summary

Here we describe the modules of Optimus; [the code](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) and [library of tasks](https://github.com/HumanCellAtlas/skylab/tree/master/library/tasks) are available through Github.

Overall, the workflow:
1. Converts R2 Fastq file to BAM
2. Corrects and attaches 10X Barcodes 
3. Aligns reads to the genome
4. Annotates genes with aligned reads
5. Generates an expression count matrix in a UMI-aware fashion
6. Detects empty droplets
7. Calculates summary statistics (metrics) and count
8. Produces a count matrix
9. Returns output in BAM, Zarr, or loom file formats

Special care is taken to flag but avoid the removal of reads that are not aligned or that do not contain recognizable barcodes. This design (which differs from many pipelines currently available) allows use of the entire dataset by those who may want to use alternative filtering or leverage the data for methodological development associated with the data processing.

### 2. Converting R2 Fastq file to BAM

Because the pipeline processing steps require a BAM file format, the first step of Optimus is to [convert](https://broadinstitute.github.io/picard/command-line-overview.html#FastqToSam) the R2 Fastq file, containing the alignable genomic information, to a BAM file.

### 2. Correcting and Attaching Cell Barcodes

Although the function of the cell barcodes is to identify unique cells, barcode errors can arise during sequencing (such as incorporation of the barcode into contaminating DNA or sequencing and PCR errors), making it difficult to distinguish unique cells from artifactual appearances of the barcode. Barcode errors are evaluated in the [Attach10xBarcodes](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/Attach10xBarcodes.wdl) step mentioned above, which compares the sequences against a whitelist of known barcode sequences.

Next, the [Attach10xBarcodes](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/Attach10xBarcodes.wdl) step appends the UMI and Cell Barcode sequences from R1 to the corresponding R2 sequence as tags, in order to properly label the genomic information for alignment.

The output file contains the reads with correct barcodes, including barcodes that came within one edit distance ([Levenshtein distance](http://www.levenshtein.net/)) of matching the whitelist of barcode sequences and were corrected by this tool. Correct barcodes are assigned a “CB” tag. Uncorrectable barcodes (with more than one error) are preserved and given a “CR” (Cell barcode Raw) tag. Cell barcode quality scores are also preserved in the file under the “CY” tag.

The various BAM files are then [scattered](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/ScatterBam.wdl) and [split](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/SplitBamByCellBarcode.wdl) into groups according to cell barcode to facilitate the following processing steps. 

### 3. Alignment

The [STAR alignment](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/StarAlignBamSingleEnd.wdl) software ([Dobin, et al., 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) is used to map barcoded reads in the BAM file to the human genome primary assembly reference (see table above for version information). STAR (Spliced Transcripts Alignment to a Reference) is widely-used for RNA-seq alignment and identifies the best matching location(s) on the reference for each sequencing read.

### 4. Gene Annotation

The [TagGeneExon](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/TagGeneExon.wdl) tool then annotates each read with the type of sequence to which it aligns. These annotations include INTERGENIC, INTRONIC, and EXONIC, and are stored using the XF BAM tag. In cases where the gene corresponds to an intron or exon, the name of the gene that overlaps the alignment is associated with the read and stored using the GE BAM tag.

### 5. UMI Correction

UMIs are designed to distinguish unique transcripts present in the cell at lysis from those arising from PCR amplification of these same transcripts. But, like cell barcodes, UMIs can also be incorrectly sequenced or amplified. Optimus uses the [UMI-tools software package](https://github.com/CGATOxford/UMI-tools), which applies a network-based method to account for such errors ([Smith, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/)). Optimus uses the “directional” method.

### 6. Identification of Empty Droplets

In addition, the pipeline runs the EmptyDrops function from the [dropletUtils](http://bioconductor.org/packages/release/bioc/html/DropletUtils.html) R package to identify cell barcodes that correspond to empty droplets. Empty droplets are those that did not encapsulate a cell but instead acquired cell-free RNA from the solution in which the cells resided -- such as secreted RNA or RNA released when some cells lysed in solution ([Lun, et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/?term=30902100)). This ambient RNA can serve as a substrate for reverse transcription, leading to a small number of background reads. Cell barcodes that are not believed to represent cells are identified in the metrics and raw information from dropletUtils is provided to the user.

### 7. Metric calculation

A number of [quality control tools](https://github.com/HumanCellAtlas/sctools) are used to assess the quality of the data output each time this pipeline is run. For a list of the tools and information about each one please see our [QC Metrics](/pipelines/hca-pipelines/data-processing-pipelines/qc-metrics) page.

### 8. Count Matrix Construction

The pipeline outputs a count matrix that contains, for each cell barcode and for each gene, the number of molecules that were observed. The script that generates this matrix evaluates every read. It discards any read that maps to more than one gene, and counts any remaining reads provided the triplet of cell barcode, molecule barcode, and gene name is unique, indicating the read originates from a single transcript present at the time of lysis of the cell.


### 9. Outputs
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
