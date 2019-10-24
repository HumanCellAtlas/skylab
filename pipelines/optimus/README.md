| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| optimus_v1.3.6 | October 10, 2019 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) | Please file GitHub issues in skylab or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |

# Table of Contents
- [Optimus Pipeline Overview](#optimus-pipeline-overview)
  * [Introduction to the Optimus Workflow](#introduction-to-the-optimus-workflow)
  * [Quick Start Table](#quick-start-table)
- [Set-up](#set-up)
  * [Optimus Installation and Requirements](#optimus-installation-and-requirements)
  * [Inputs](#inputs)
    + [Sample Data Input](#sample-data-input)
    + [Additional Reference Inputs](#additional-reference-inputs)
- [Running Optimus](#running-optimus)
  * [Optimus Modules Summary](#optimus-modules-summary)
    + [1. Converting R2 Fastq File to UBAM](#1-converting-r2-fastq-file-to-ubam)
    + [2. Correcting and Attaching Cell Barcodes](#2-correcting-and-attaching-cell-barcodes)
    + [3. Alignment](#3-alignment)
    + [4. Gene Annotation](#4-gene-annotation)
    + [5. UMI Correction](#5-umi-correction)
    + [6. Summary Metric Calculation](#6-summary-metric-calculation)
    + [7. Expression Matrix Construction](#7-expression-matrix-construction)
    + [8. Identification of Empty Droplets](#8-identification-of-empty-droplets)
    + [9. Outputs](#9-outputs)
- [Versioning](#versioning)
- [Additional Notes](#additional-notes)

# Optimus Pipeline Overview
![Diagram](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus_diagram.png)

## Introduction to the Optimus Workflow

Optimus is a pipeline developed by the Data Coordination Platform (DCP) of the [Human Cell Atlas (HCA) Project](https://data.humancellatlas.org/) that supports processing of any 3' single-cell expression data generated with the [10X Genomic V2 and V3 assay](https://www.10xgenomics.com/solutions/single-cell/). It is an alignment and transcriptome quantification pipeline that corrects Cell Barcodes, aligns reads to the genome, corrects Unique Molecular Identifiers (UMIs), generates an expression matrix in a UMI-aware manner, detects empty droplets, calculates summary metrics for genes and cells, returns read outputs in BAM format, and returns cell gene expression in numpy matrix, Zarr, and Loom file formats. Special care is taken to keep all reads that may be useful to the downstream user, such as unaligned reads or reads with uncorrectable barcodes. This design provides flexibility to the downstream user and allows for alternative filtering or leveraging the data for novel methodological development.

Optimus has been validated for analyzing both [human](https://github.com/HumanCellAtlas/skylab/blob/master/benchmarking/optimus/optimus_report.rst) and [mouse](https://docs.google.com/document/d/1_3oO0ZQSrwEoe6D3GgKdSmAQ9qkzH_7wrE7x6_deL10/edit) data sets. More details about the human validation can be found in the [in the original file](https://docs.google.com/document/d/158ba_xQM9AYyu8VcLWsIvSoEYps6PQhgddTr9H0BFmY/edit).

## Quick Start Table

| Pipeline Features | Description | Source |
|-------------------|---------------------------------------------------------------|-----------------------|
| Assay Type | 10x Single Cell Expression (v2 and v3) | [10x Genomics](https://www.10xgenomics.com)
| Overall Workflow  | Quality control module and transcriptome quantification module | Code available from [Github](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) |
| Workflow Language | WDL          | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence|GRCh38 human genome primary sequence and M21 (GRCm38.p6) mouse genome primary sequence|GENCODE [Human](https://www.gencodegenes.org/human/release_27.html) and [Mouse](https://www.gencodegenes.org/mouse/release_M21.html) 
| Transcriptomic Reference Annotation | V27 GenCode human transcriptome and M21 mouse transcriptome | GENCODE [Human](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz) and [Mouse](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz) |
| Aligner           | STAR (v.2.5.3)       | [Dobin, et al.,2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) |
| Transcript Quantification | Utilities for processing large-scale single cell datasets | [Sctools](https://github.com/HumanCellAtlas/sctools)                          
| Data Input File Format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |                     
| Data Output File Format | File formats in which Optimus output is provided | [BAM](http://samtools.github.io/hts-specs/), [Zarr version 2](https://zarr.readthedocs.io/en/stable/spec/v2.html), Python numpy arrays (internal), Loom (generated with [Loompy v.2.0.17)](http://loompy.org/) |

# Set-up

## Optimus Installation and Requirements
The Optimus pipeline code can be downloaded by cloning the GitHub repository [Skylab](https://github.com/HumanCellAtlas/skylab/). For the latest release of Optimus, please see the realease tags prefixed with "optimus" [here](https://github.com/HumanCellAtlas/skylab/releases). 

Optimus can be deployed using [Cromwell](https://software.broadinstitute.org/wdl/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. Optimus can also be run in [Terra](https://app.terra.bio/#workspaces/help-gatk/HCA_Optimus_Pipeline), a cloud-based analysis platform. In this featured workspace the user will find the Optimus pipeline, configurations, required reference data and other inputs, and example testing data.

## Inputs

Optimus pipeline inputs are detailed in a json file, such as in this [example](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/example_test_inputs.json). 

### Sample Data Input

Each 10X v2 and v3 3’ sequencing experiment generates triplets of fastq files for any given sample:  

1. A forward reads (r1_fastq), containing the unique molecular identifier (UMI) and cell barcode sequences
2. A reverse reads (r2_fastq), which contain the alignable genomic information from the mRNA transcript 
3. An index fastq (i1_fastq) that contains the sample barcodes, when provided by the sequencing facility

Note: Optimus is currently a single sample pipeline, but can take in multiple sets of fastqs for a sample that has been split over lanes of sequencing. 

### Additional Reference Inputs

The json file also contains metadata for the following reference information:

* Whitelist: a list of known cell barcodes from [10X genomics](https://www.10xgenomics.com/)
* Tar_star_reference: TAR file containing a specifes-specific reference genome and gtf; it is generated using the [StarMkRef.wdl](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/StarMkref.wdl)
* Sample_id: a unique name describing the biological sample or replicate that corresponds with the original fastq files
* Annotations_gtf: a GTF containing gene annotations used for gene tagging (must match gft in STAR reference)
* Chemistry: an optional description of whether data was generated with V2 or V3 chemistry


# Running Optimus

* The [Optimus.wdl](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) in the pipelines/optimus folder of the HCA Skylab repository implements the workflow by importing individual modules ("tasks" written in  WDL) from the Skylab [Library](https://github.com/HumanCellAtlas/skylab/tree/master/library) folder.


## Optimus Modules Summary

Here we describe the modules ("tasks") of the Optimus Pipeliene; [the code](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) and [library of tasks](https://github.com/HumanCellAtlas/skylab/tree/master/library/tasks) are available through GitHub.

Overall, the workflow:
1. Converts R2 fastq file (containing alignable genomic information) to an unaligend BAM (UBAM)
2. Corrects and attaches 10X Barcodes using the R1 Fastq file 
3. Aligns reads to the genome with STAR v.2.5.3a
4. Annotates genes with aligned reads
5. Corrects UMIs
6. Calculates summary metrics
7. Produces a UMI-aware expression matrix
8. Detects empty droplets
9. Returns a GA4GH compliant BAM and metric matrix in Zarr or Loom formats

Special care is taken to flag but avoid the removal of reads that are not aligned or that do not contain recognizable barcodes. This design (which differs from many pipelines currently available) allows use of the entire dataset by those who may want to use alternative filtering or leverage the data for methodological development associated with the data processing.

### 1. Converting R2 Fastq File to UBAM

Unlike fastq files, BAM files enable researchers to keep track of important metadata throughout all data processing steps. The first step of Optimus is to [convert](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/FastqToUBam.wdl) the R2 fastq file, containing the alignable genomic information, to an unaligned BAM (UBAM) file.

### 2. Correcting and Attaching Cell Barcodes

Although the function of the cell barcodes is to identify unique cells, barcode errors can arise during sequencing (such as incorporation of the barcode into contaminating DNA or sequencing and PCR errors), making it difficult to distinguish unique cells from artifactual appearances of the barcode. The [Attach10xBarcodes](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/Attach10xBarcodes.wdl) task uses [sc-tools v.0.3.4](https://github.com/HumanCellAtlas/sctools) to evaluate barcode errors by comparing the R1 fastq sequences against a whitelist of known barcode sequences. The task then appends the UMI and Cell Barcode sequences from the R1 fastq to the UBAM sequence as tags, properly labeling the genomic information for alignment. 

The output is a UBAM file containing the reads with correct barcodes, including barcodes that came within one edit distance ([Levenshtein distance](http://www.levenshtein.net/)) of matching the whitelist of barcode sequences and were corrected by this tool. Correct barcodes are assigned a “CB” tag. Uncorrectable barcodes (with more than one error) are preserved and given a “CR” (Cell barcode Raw) tag. Cell barcode quality scores are also preserved in the file under the “CY” tag.

To facilitate subsequent processing steps, the pipeline then [scatters](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/ScatterBam.wdl) and [splits](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/SplitBamByCellBarcode.wdl) the corrected UBAM files into groups according to cell barcode. 

### 3. Alignment

Optimus uses the [STAR alignment](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/StarAlignBamSingleEnd.wdl) task to map barcoded reads in the UBAM file to the genome primary assembly reference (see table above for version information). This task uses STAR (Spliced Transcripts Alignment to a Reference; [Dobin, et al., 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) a standard, splice-aware, RNA-seq alignment tool.

### 4. Gene Annotation

The [TagGeneExon](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/TagGeneExon.wdl) task then uses [Drop-seq tools v.1.13](https://github.com/broadinstitute/Drop-seq) to annotate each read with the type of sequence to which it aligns. These annotations include INTERGENIC, INTRONIC, and EXONIC, and are stored using the XF BAM tag. In cases where the gene corresponds to an intron or exon, the name of the gene that overlaps the alignment is associated with the read and stored using the GE BAM tag.


### 5. UMI Correction

UMIs are designed to distinguish unique transcripts present in the cell at lysis from those arising from PCR amplification of these same transcripts. But, like cell barcodes, UMIs can also be incorrectly sequenced or amplified. Optimus uses the [UmiCorrection task](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/UmiCorrection.wdl) to apply a network-based, "directional" method ([Smith, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/)) to account for such errors using [Umi-tools v.0.0.1](https://pypi.org/project/umi-tools/0.0.1/) 

### 6. Summary Metric Calculation

The Metrics task calls the [SequenceDataWithMoleculeTagMetrics.wdl](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/SequenceDataWithMoleculeTagMetrics.wdl) to calculate summary metrics which are used assess the quality of the data output each time this pipeline is run. This task uses sctools v.0.3.3](https://github.com/HumanCellAtlas/sctools). These metrics are included in the ZARR and [Loom](link to Loom schema) output files.

### 7. Expression Matrix Construction

The Optimus [Count]((https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/CreateCountMatrix.wdl)) task evaluates every read in the BAM file and creates a UMI-aware expression matrix using [Drop-seq tools](https://github.com/broadinstitute/Drop-seq). This matrix contains the number of molecules that were observed for each cell barcode and for each gene. The task discards any read that maps to more than one gene, and counts any remaining reads provided the triplet of cell barcode, molecule barcode, and gene name is unique, indicating the read originates from a single transcript present at the time of cell lysis.

### 8. Identification of Empty Droplets

Empty droplets are lipid droplets that did not encapsulate a cell during 10X sequencing, but instead acquired cell-free RNA (secreted RNA or RNA released during cell lysis ([Lun, et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/?term=30902100)) from the solution in which the cells resided. This ambient RNA can serve as a substrate for reverse transcription, leading to a small number of background reads. The Optimus pipeline calls the RunEmptyDrops task which uses the [dropletUtils v.0.1.1](http://bioconductor.org/packages/release/bioc/html/DropletUtils.html) R package to flag cell barcodes that represent empty droplets rather than cells. These metrics are stored in the output Zarr and [Loom](link to Loom schema) files. 

### 9. Outputs

Output files of the pipeline include:

1. Cell x Gene unnormalized expression matrix
2. Unfiltered, sorted BAM file with [barcode and downstream analysis Tags](link to UBAM Tag Description)
3. Cell metadata, including cell metrics
4. Gene metadata, including gene metrics

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

# Versioning

| Optimus Release Version | Date | Release Note | 
| :---: | :---: | :---: |
| v1.3.6 (current) | 09/23/2019 | Optimus now optionally outputs a Loom formatted count matrix, with the default being true. |
| v1.3.3 | 08/29/2019 | This version and newer have been validated to additionally support Mouse data. The gene expression per cell is now counted by gencode geneID instead of gene name. There is an additional output mapping geneID to gene name provided. This is a breaking change. | 
| v1.0.0 |03/30/2019 | Initial pipeline release. Validated on hg38 gencodev27. | 

# Have Suggestions? 

Coming soon, we will have a document dedicated to open issues! Please help us make our tools better. In the mean time, please contact Kylee Degatano](mailto:kdegatano@broadinstitute.org) for pipeline-related questions.




