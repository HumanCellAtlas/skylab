# Run SmartSeq2 scRNA-Seq pipeline on Batch of Samples
This workflow is for running SmartSeq RNA-Seq pipeline through batch of scRNA-Seq dataset.
## Features
 - Aligns reads of each sample to reference genome (using STAR)
 - Performs quality control on generated BAM files (using Picard)
   * [RNA metrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics)
   * [Alignment metrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html#AlignmentSummaryMetrics)
   * [Duplication metrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html#DuplicationMetrics)
   * [Insertion Size metrics](https://broadinstitute.github.io/picard/picard-metric-definitions.html#InsertSizeMetrics)
   * [RSEM transcrptome quantification](http://deweylab.biostat.wisc.edu/rsem/rsem-calculate-expression.html)
   * [FeatureCounts gene/exon/transcript quantification](http://bioinf.wehi.edu.au/featureCounts/). FeatureCounts counts both unique mapping reads and multiple mapping reads.  

# Getting Started
## Download
Use git clone `git clone git@github.com:HumanCellAtlas/skylab.git`
## Requirements
### Dockers
- STAR(2_5_3a) docker: `humancellatlas/star_dev:v1`
- RSEM docker: `humancellatlas/rsem`
- Picard docker: `broadinstitute/genomes-in-the-cloud:2.3.1-1504795437`
### Reference Genome And Gene Annotaiton
 - Star reference indexes `star.tar`
 - RSEM reference indexes `rsem.tar`
 - Genome fasta
 - refFlat file
 - rRNA interval list
 - Gencode annotation(gtf)
## Input Json
Main workflow `ss2_multi_sample_wf.wdl` takes a json file as input. An example json file is listed below:
```json
{
  "Ss2RunMultiSample.ref_fasta": "reference fasta File",
  "Ss2RunMultiSample.rsem_genome": "tarball file of rsem reference indexes",
  "Ss2RunMultiSample.ref_flat": "refflat file",
  "Ss2RunMultiSample.sra_dir": "the location to store data extract from sra/geo",
  "Ss2RunMultiSample.star_genome": "tarball file of star reference indexes",
  "Ss2RunMultiSample.sra_list_file": "list of sra ids",
  "Ss2RunMultiSample.rrna_interval": "Ribosomal coordination in  interval_list format",
  "Ss2RunMultiSample.gtf": "a gtf file of gene annotaiton"
}
```
### Dependencies
Main workflow imports a subworkflow `ss2_single_sample_wf.wdl`
## Test and Demo
A small demo dataset(see in `ss2_multi_sample_wf_demo.json`) can be used to test this workflow. Reference used in this demo dataset is built from chr21 only.

