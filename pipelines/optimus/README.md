| Pipeline Version | Date Updated | Documentation Author | Questions or Feedback |
| :----: | :---: | :----: | :--------------: |
| [optimus_v1.4.0](https://github.com/HumanCellAtlas/skylab/releases/tag/optimus_v1.4.0) | November 08, 2019 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) | Please file GitHub issues in skylab or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |
| [optimus_v1.3.6](https://github.com/HumanCellAtlas/skylab/releases/tag/optimus_v1.3.6) | October 10, 2019 | [Elizabeth Kiernan](mailto:ekiernan@broadinstitute.org) | Please file GitHub issues in skylab or contact [Kylee Degatano](mailto:kdegatano@broadinstitute.org) |

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
- [Pipeline Improvements](#have-suggestions)

# Optimus Pipeline Overview
![Diagram](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus_diagram.png)

## Introduction to the Optimus Workflow

Optimus is a pipeline developed by the Data Coordination Platform (DCP) of the [Human Cell Atlas (HCA) Project](https://data.humancellatlas.org/) that supports processing of any 3' single-cell expression data generated with the [10x Genomic V2 or V3 assay](https://www.10xgenomics.com/solutions/single-cell/). It is an alignment and transcriptome quantification pipeline that corrects Cell Barcodes, aligns reads to the genome, corrects Unique Molecular Identifiers (UMIs), generates an expression matrix in a UMI-aware manner, calculates summary metrics for genes and cells, detects empty droplets, returns read outputs in BAM format, and returns cell gene expression in numpy matrix, Zarr, and Loom file formats. Special care is taken to keep all reads that may be useful to the downstream user, such as unaligned reads or reads with uncorrectable barcodes. This design provides flexibility to the downstream user and allows for alternative filtering or leveraging the data for novel methodological development.

Optimus has been validated for analyzing both [human](https://github.com/HumanCellAtlas/skylab/blob/master/benchmarking/optimus/optimus_report.rst) and [mouse](https://docs.google.com/document/d/1_3oO0ZQSrwEoe6D3GgKdSmAQ9qkzH_7wrE7x6_deL10/edit) data sets. More details about the human validation can be found in the [in the original file](https://docs.google.com/document/d/158ba_xQM9AYyu8VcLWsIvSoEYps6PQhgddTr9H0BFmY/edit).

## Quick Start Table

| Pipeline Features | Description | Source |
|-------------------|---------------------------------------------------------------|-----------------------|
| Assay Type | 10x Single Cell Expression (v2 and v3) | [10x Genomics](https://www.10xgenomics.com)
| Overall Workflow  | Quality control module and transcriptome quantification module | Code available from [Github](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) |
| Workflow Language | WDL | [openWDL](https://github.com/openwdl/wdl) |
| Genomic Reference Sequence| GRCh38 human genome primary sequence and M21 (GRCm38.p6) mouse genome primary sequence | GENCODE [Human](https://www.gencodegenes.org/human/release_27.html) and [Mouse](https://www.gencodegenes.org/mouse/release_M21.html) 
| Transcriptomic Reference Annotation | V27 GENCODE human transcriptome and M21 mouse transcriptome | GENCODE [Human](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz) and [Mouse](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gff3.gz) |
| Aligner  | STAR (v.2.5.3a) | [Dobin, et al.,2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) |
| Transcript Quantification | Utilities for processing large-scale single cell datasets | [Sctools](https://github.com/HumanCellAtlas/sctools)                          
| Data Input File Format | File format in which sequencing data is provided | [FASTQ](https://academic.oup.com/nar/article/38/6/1767/3112533) |                     
| Data Output File Format | File formats in which Optimus output is provided | [BAM](http://samtools.github.io/hts-specs/), [Zarr version 2](https://zarr.readthedocs.io/en/stable/spec/v2.html), Python numpy arrays (internal), Loom (generated with [Loompy v.2.0.17)](http://loompy.org/) |

# Set-up

## Optimus Installation and Requirements
The Optimus pipeline code can be downloaded by cloning the GitHub repository [skylab](https://github.com/HumanCellAtlas/skylab/). For the latest release of Optimus, please see the release tags prefixed with "optimus" [here](https://github.com/HumanCellAtlas/skylab/releases). 

Optimus can be deployed using [Cromwell](https://software.broadinstitute.org/wdl/), a GA4GH compliant, flexible workflow management system that supports multiple computing platforms. Optimus can also be run in [Terra](https://app.terra.bio/#workspaces/help-gatk/HCA_Optimus_Pipeline), a cloud-based analysis platform. In this featured workspace, the user will find the Optimus pipeline, configurations, required reference data and other inputs, and example testing data.

## Inputs

Optimus pipeline inputs are detailed in a json file, such as in this [example](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/example_test_inputs.json). 

### Sample Data Input

Each 10x v2 and v3 3’ sequencing experiment generates triplets of fastq files for any given sample:  

1. A forward reads (r1_fastq), containing the unique molecular identifier (UMI) and cell barcode sequences
2. A reverse reads (r2_fastq), which contain the alignable genomic information from the mRNA transcript 
3. An index fastq (i1_fastq) that contains the sample barcodes, when provided by the sequencing facility

Note: Optimus is currently a single sample pipeline, but can take in multiple sets of fastqs for a sample that has been split over multiple lanes of sequencing. Additionally, Optimus does not support demultiplexing even though it accepts index fastq files. 

### Additional Reference Inputs

The json file also contains metadata for the following reference information:

* Whitelist: a list of known cell barcodes from [10x genomics](https://www.10xgenomics.com/) that corresponds to the V2 or V3 chemistry
* Tar_star_reference: TAR file containing a species-specific reference genome and gtf; it is generated using the [StarMkRef.wdl](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/StarMkref.wdl)
* Sample_id: a unique name describing the biological sample or replicate that corresponds with the original fastq files
* Annotations_gtf: a GTF containing gene annotations used for gene tagging (must match GTF in STAR reference)
* Chemistry: an optional string description of whether data was generated with 10x V2 or V3 chemistry
  * Optional string: "tenX_v2" (default) or "tenX_v3"
   * Note: Optimus validates this string. If the string does not match these options, the pipeline will fail. You can remove the checks by setting "force_no_check = true" in the input json.
# Running Optimus

* The [Optimus.wdl](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) in the pipelines/optimus folder of the HCA skylab repository implements the workflow by importing individual modules ("tasks" written in  WDL script) from the skylab [Library](https://github.com/HumanCellAtlas/skylab/tree/master/library) folder.

## Optimus Modules Summary

Here we describe the modules ("tasks") of the Optimus pipeline; [the code](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) and [library of tasks](https://github.com/HumanCellAtlas/skylab/tree/master/library/tasks) are available through GitHub.

Overall, the workflow:
1. Converts R2 fastq file (containing alignable genomic information) to an unaligned BAM (UBAM)
2. Corrects and attaches 10x Barcodes using the R1 Fastq file 
3. Aligns reads to the genome with STAR v.2.5.3a
4. Annotates genes with aligned reads
5. Corrects UMIs
6. Calculates summary metrics
7. Produces a UMI-aware expression matrix
8. Detects empty droplets
9. Returns a GA4GH compliant BAM and an expression matrix in Zarr or Loom formats

The tools each Optimus task employs are detailed in the following table:

| Task | Tool | 
| --- | --- |
| [Attach10xBarcodes](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/Attach10xBarcodes.wdl) |	[sctools](https://sctools.readthedocs.io/en/latest/sctools.html) |
| [StarAlignBamSingleEnd](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/StarAlignBamSingleEnd.wdl) |	[STAR](https://github.com/alexdobin/STAR) |
| [TagGeneExon](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/TagGeneExon.wdl) |	[Drop-seq](https://github.com/broadinstitute/Drop-seq) |
| [UmiCorrection](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/UmiCorrection.wdl) |	[Umi-tools](https://github.com/CGATOxford/UMI-tools) |
| [SequenceDataWithMoleculeTagMetrics](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/SequenceDataWithMoleculeTagMetrics.wdl) |	[sctools](https://sctools.readthedocs.io/en/latest/sctools.html) |
| [RunEmptyDrops](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/RunEmptyDrops.wdl) |	[dropletUtils](https://bioconductor.org/packages/release/bioc/html/DropletUtils.html) |
| [CreateCountMatrix](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/CreateCountMatrix.wdl) |	[Drop-seq](https://github.com/broadinstitute/Drop-seq) |
| [FastqToUBAM](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/FastqToUBam.wdl)	| [picard](https://github.com/broadinstitute/picard) |
| [SplitBamByCellBarcode](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/SplitBamByCellBarcode.wdl) |	[sctools](https://sctools.readthedocs.io/en/latest/sctools.html) |
| [TagSortBam](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/TagSortBam.wdl) |	[sctools](https://sctools.readthedocs.io/en/latest/sctools.html) |
| [Picard](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/Picard.wdl)	| [picard](https://github.com/broadinstitute/picard) |

The Optimus pipeline takes special care to flag but avoid the removal of reads that are not aligned or that do not contain recognizable barcodes. This design (which differs from many pipelines currently available) allows the use of the entire dataset by those who may want to use alternative filtering or leverage the data for methodological development associated with the data processing. More information about the different tags used to flag the data can be found [here](Bam_tags.md).

### 1. Converting R2 Fastq File to UBAM

Unlike fastq files, BAM files enable researchers to keep track of important metadata throughout all data processing steps. The first step of Optimus is to [convert](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/FastqToUBam.wdl) the R2 fastq file, containing the alignable genomic information, to an unaligned BAM (UBAM) file.

### 2. Correcting and Attaching Cell Barcodes

Although the function of the cell barcodes is to identify unique cells, barcode errors can arise during sequencing (such as incorporation of the barcode into contaminating DNA or sequencing and PCR errors), making it difficult to distinguish unique cells from artifactual appearances of the barcode. The [Attach10xBarcodes](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/Attach10xBarcodes.wdl) task uses [sctools](https://github.com/HumanCellAtlas/sctools) to evaluate barcode errors by comparing the R1 fastq sequences against a whitelist of known barcode sequences. The task then appends the UMI and Cell Barcode sequences from the R1 fastq to the UBAM sequence as [tags](Bam_tags.md). 

The output is a UBAM file containing the reads with corrected barcodes, including barcodes that came within one edit distance ([Levenshtein distance](http://www.levenshtein.net/)) of matching the whitelist of barcode sequences and were corrected by this tool. Correct barcodes are assigned a “CB” tag. Uncorrectable barcodes (with more than one error) are preserved and given a “CR” (Cell barcode Raw) tag. Cell barcode quality scores are also preserved in the file under the “CY” tag.

To enable parallelization, the pipeline then [scatters](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/ScatterBam.wdl) and [splits](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/SplitBamByCellBarcode.wdl) the corrected UBAM files into groups according to cell barcode. 

### 3. Alignment

Optimus uses the [STAR alignment](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/StarAlignBamSingleEnd.wdl) task to map barcoded reads in the UBAM file to the genome primary assembly reference (see table above for version information). This task uses STAR (Spliced Transcripts Alignment to a Reference) a standard, splice-aware, RNA-seq alignment tool [(Dobin, et al., 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/).

### 4. Gene Annotation

The [TagGeneExon](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/TagGeneExon.wdl) task then uses [Drop-seq tools](https://github.com/broadinstitute/Drop-seq) to [annotate each read](https://github.com/HumanCellAtlas/skylab/blob/LK_BAM_TAGS/pipelines/optimus/Optimus_BAM_tags.md) with the type of sequence to which it aligns. These annotations include INTERGENIC, INTRONIC, and EXONIC, and are stored using the XF BAM [tag](Bam_tags.md). In cases where the gene corresponds to an intron or exon, the name of the gene that overlaps the alignment is associated with the read and stored using the GE BAM tag.

### 5. UMI Correction

UMIs are designed to distinguish unique transcripts present in the cell at lysis from those arising from PCR amplification of these same transcripts. But, like cell barcodes, UMIs can also be incorrectly sequenced or amplified. The [UmiCorrection] (https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/UmiCorrection.wdl)task uses [Umi-tools v.0.0.1](https://pypi.org/project/umi-tools/0.0.1/) to apply a network-based, "directional" correction method ([Smith, et al., 2017](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5340976/)) to account for such errors. 

### 6. Summary Metric Calculation

The [Metrics](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/SequenceDataWithMoleculeTagMetrics.wdl) task uses [sctools](https://github.com/HumanCellAtlas/sctools) to calculate summary metrics which help assess the quality of the data output each time this pipeline is run. These metrics are included in the Zarr and [Loom](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Loom_schema.md) output files.

### 7. Expression Matrix Construction

The Optimus [Count](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/CreateCountMatrix.wdl) task evaluates every read in the BAM file and creates a UMI-aware expression matrix using [Drop-seq tools](https://github.com/broadinstitute/Drop-seq). This matrix contains the number of molecules that were observed for each cell barcode and for each gene. The task discards any read that maps to more than one gene, and counts any remaining reads provided the triplet of cell barcode, molecule barcode, and gene name is unique, indicating the read originates from a single transcript present at the time of cell lysis.

### 8. Identification of Empty Droplets

Empty droplets are lipid droplets that did not encapsulate a cell during 10x sequencing, but instead acquired cell-free RNA (secreted RNA or RNA released during cell lysis) from the solution in which the cells resided ([Lun, et al., 2018](https://www.ncbi.nlm.nih.gov/pubmed/?term=30902100). This ambient RNA can serve as a substrate for reverse transcription, leading to a small number of background reads. The Optimus pipeline calls the [RunEmptyDrops](https://github.com/HumanCellAtlas/skylab/blob/master/library/tasks/RunEmptyDrops.wdl) task which uses the [dropletUtils v.0.1.1](http://bioconductor.org/packages/release/bioc/html/DropletUtils.html) R package to flag cell barcodes that represent empty droplets rather than cells. A cell will be flagged if it contains fewer than 100 molecules. These metrics are stored in the output Zarr and [Loom](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Loom_schema.md) files. 

### 9. Outputs

Output files of the pipeline include:

1. Cell x Gene unnormalized, but UMI-corrected, expression matrices
2. Unfiltered, sorted BAM file with barcode and downstream analysis [tags](Bam_tags.md)
3. Cell metadata, including cell metrics
4. Gene metadata, including gene metrics

Following are the types of files produced from the pipeline.

| Output Name | Filename, if applicable | Output Type |Output Format |
| ------ |------ | ------ | ------ | 
| pipeline_version | | Version of the processing pipeline run on this data | String | 
| bam | merged.bam | aligned bam | bam |
| matrix_row_index | sparse_counts_row_index.npy | Index of cells in expression matrix | Numpy array index |
| matrix_col_index | sparse_counts_col_index.npy | Index of genes in expression matrix | Numpy array index | 
| cell_metrics | merged-cell-metrics.csv.gz | cell metrics | compressed csv | Matrix of metrics by cells | 
| gene_metrics | merged-gene-metrics.csv.gz | gene metrics | compressed csv | Matrix of metrics by genes | 
| zarr_output_files | {unique_id}.zarr!.zattrs | Zarr | Array | 
| loom_output_file | output.loom | Loom | Loom | Loom file with expression data and metadata | N/A |


The Zarr array is the default output. The Zarr schema version is detailed in the array as 'optimus_output_schema_version'. The schema version is specified to the Zarr using the [create_zarr_optimus.py](https://github.com/HumanCellAtlas/skylab/blob/master/docker/zarr-output/create_zarr_optimus.py) script. 

The Loom file is an optional output that is specified in the "meta" section of the [Optimus workflow](https://github.com/HumanCellAtlas/skylab/blob/master/pipelines/optimus/Optimus.wdl) with the following boolean command:

> "Boolean output_loom = false"

To obtain a Loom file, the boolean parameter "false" must be changed to "true". 

# Versioning

| Optimus Release Version | Date | Release Note | 
| :---: | :---: | :---: |
| [v1.4.0](https://github.com/HumanCellAtlas/skylab/releases/tag/optimus_v1.4.0) (current) | 11/08/2019 | Optimus now supports processing 10X V3 Chemistry scRNAseq data. By default, the pipeline processes 10X V2, V3 must be specified as an input. EmptyDrops calls are now a column in Zarr output. |
| v1.3.6 | 09/23/2019 | Optimus now optionally outputs a Loom formatted count matrix. |
| v1.3.3 | 08/29/2019 | This version and newer have been validated to additionally support Mouse data. The gene expression per cell is now counted by gencode geneID instead of gene name. There is an additional output mapping geneID to gene name provided. This is a breaking change. | 
| v1.0.0 |03/30/2019 | Initial pipeline release. Validated on hg38 GENCODE v27. | 

# Have Suggestions? 

Coming soon, we will have a GitHub document dedicated to open issues! In the meantime, please help us make our tools better by contacting [Kylee Degatano](mailto:kdegatano@broadinstitute.org) for pipeline-related suggestions or questions.


