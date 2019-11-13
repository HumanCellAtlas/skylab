# Snap-ATAC pipeline Overview

The snap-atac pipeline is a pipeline for processing scATAC-seq datasets, based on the pipeline that is available (in the snap-atac repository)[https://github.com/r3fang/SnapATAC) and [the snap-utils repository](https://github.com/r3fang/SnapATAC).

## Pipeline Inputs

The pipeline accepts paired reads in the form of FASTQ files. The current version of the pipeline requires that the cellular barcodes that in most protocols are present in index read 1 and index read2 have already been appended to the FASTQ read file names. The following is an example of the format of the expected input. The full cell barcode must form the first part of the read name (for both R1 and R2 files) and be separated from the rest of the line by a colon.

```
@CAGTTGCACGTATAGAACAAGGATAGGATAAC:7001113:915:HJ535BCX2:1:1106:1139:1926 1:N:0:0
ACCCTCCGTGTGCCAGGAGATACCATGAATATGCCATAGAACCTGTCTCT
+
DDDDDIIIIIIIIIIIIIIHHIIIIIIIIIIIIIIIIIIIIIIIIIIIII
```

The following table outlines the inputs to the pipeline in more detail.

| input name       | input type    | description                                                                |
|------------------|---------------|----------------------------------------------------------------------------|
| input_fastq1     | File          | fastq file of the first reads (R1)                                         |
| input_fastq2     | File          | fastq file of the second reads (R2)                                        |
| genome_name      | String        | name of the genomic reference used                                         |
| input_reference  | File          | reference bundle that is generated with bwa-mk-index-wdl                   |
| output_bam       | String        | name of the output bam to generate                                         |


`input_reference` is an input reference bundle that has been generated with the bwa-mk-index-wdl accessory pipeline. This accessory pipeline can be found [here](https://github.com/HumanCellAtlas/skylab/tree/master/library/accessory_workflows/build_bwa_reference).

## Summary of pipeline execution steps

The pipeline is composed of five steps:

| step name        | step description                                                                         |
|------------------|------------------------------------------------------------------------------------------|
| AlignPairEnd     | Align the fastq files to the genome                                                      |
| SnapPre          | Initial generation of snap file                                                          |
| SnapCellByBin    | Binning of data by genomic bins                                                          |
| MakeCompliantBAM | Generation of a GA4GH compliant BAM                                                      |
| BreakoutSnap     | Extraction of tables from snap file into text format (for testing and user availability) |

## Pipeline outputs

The pipeline outputs the following files:

| output file name              | description            |
|-------------------------------|------------------------|
| output_snap_qc                | Quality control file corresponding to the snap file |
| output_snap                   | Output snap file (in hdf5 container format) |
| output_aligned_bam            | Output BAM file, compliant with GA4GH |
| breakout_barcodes             | Text file containing the 'Fragments session' barcodeLen, barcodePos fields           |
| breakout_fragments            | Text file containing the 'Fragments session' fragChrom, fragLen and fragStart fields |
| breakout_binCoordinates       | Text file with the AM section ('Cell x bin accesibility' matrix), binChrom and binStart fields |
| breakout_binCounts            | Text file with the AM section ('Cell x bin accesibility' matrix), idx, idy and count fields |
| breakout_barcodesSection      | Text file with the data from the BD section ('Barcode session' table) |

The output bins for summarization of the data are of size 10kb in the snap file.

The format of the snap file, as well as what the different section contains, is described in more detail [here (external link)](https://github.com/r3fang/SnapTools) and [here](https://github.com/r3fang/SnapTools/blob/master/docs/snap_format.docx).

The output BAM file does not contain cellular barcodes in the read names. For GA4GH compliance the barcodes are placed in the CB bam file tag.