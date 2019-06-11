# Snap-ATAC pipeline Overview

The snap-atac pipeline is a pipeline for processing scATAC-seq datasets, based on the pipeline that is available (in the snap-atac repository)[https://github.com/r3fang/SnapATAC] and (the snap-utils repository)[https://github.com/r3fang/SnapATAC].

## Pipeline Inputs

| input name       | input type    | description                                                                |
|------------------|---------------|----------------------------------------------------------------------------|
| input_fastq1     | File          | fastq file of the first reads (R1)                                         |
| input_fastq2     | File          | fastq file of the second reads (R2)                                        |
| genome_name      | String        | name of the genomic reference used, currently this can be mm10 or hg38     |
| input_reference  | File          | reference bundle that is generated with bwa-mk-index-wdl                   |
| output_bam       | String        | name of the output bam to generate                                         |


The pipeline accepts paired reads in the form of FASTQ files. The current version of the pipeline requires that the cellular barcodes that in most protocols are present in index read 1 and index read2 have already been appended to the FASTQ read file names. The following is an example of the format of the expected input. The full cell barcode must form the first part of the read name (for both R1 and R2 files) and be separated from the rest of the line by a colon.

```
@CAGTTGCACGTATAGAACAAGGATAGGATAAC:7001113:915:HJ535BCX2:1:1106:1139:1926 1:N:0:0
ACCCTCCGTGTGCCAGGAGATACCATGAATATGCCATAGAACCTGTCTCT
+
DDDDDIIIIIIIIIIIIIIHHIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

## Pipeline executation Steps
The pipeline is composed of five steps:

* AlignPairEnd
* SnapPre
* SnapCellByBin
* MakeCompliantBAM
* BreakoutSnap


## Pipeline outputs

The 