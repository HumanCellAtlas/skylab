# generate_demo_data WDL 
This pipeline extracts a subset of the reads from a full 10x run that align to a selected 
chromosome. These truncated fastq files are used as input for testing functions.
- Task 1 StarAlignSubset (using STAR program): aligns the first 100,000 reads from the second 
  illumina fastq read, which contains only genomic information (and no barcodes). 
- Task 2 (using scsequtil python library): identifies the indices of the first 10,000 reads that 
  align to the desired chromosome, as well as 2000 reads that do not align.  
- Task 3 (using scsequtil python library): extracts the reads identified in Task 2 from the 
  original fastq files, producing a truncated output usable for testing. 

## Download
Use git clone: 
 
```
git clone git@github.com:HumanCellAtlas/skylab.git
cd skylab/10x/generate_reference_bundle
```

# Requirements
## Dockers
- STAR v2_5_3a docker: `quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-1.0.0`
- python 3 scientific docker: `quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.5`

## File Inputs
- `fastq_r1`: fastq read 1 for a 10x experiment
- `fastq_r2`: fastq read 2 for the matching experiment
- `fastq_i1`: fastq index read 1 for the matching 10x experiment
- `star_genome`: STAR index for the organism matching the 10x experiment. These indices can be
  created with `BuildStarReferenceBundle.wdl` in `skylab`
- `gtf`: The gtf annotation file used to construct the `star_genome` provided above. 
- `chromosome`: the chromosome `[1-22, X, Y, M]` for reads to be extracted from.

## Runtime requirements
- Memory: 8GB
- Processors: 1
- Disk Space: 10 GB (including docker images)
- Expected time: ~ 30s

# Example Input Data
Example input data is derived from a 10 million read trunction of the 10x public peripheral 
blood mononuclear cell 
<a href=https://support.10xgenomics.com/single-cell-gene-expression/datasets>8000 cell dataset</a>. 
The remaining parameters extract a demo dataset from chromosome 21. The 
genome provided was created with`BuildStarReferenceDemo.wdl` and derived from a
human genome and annotation that was downloaded from 
<a href=https://www.gencodegenes.org/releases/current.html>GENCODE</a>

# Output Description
- `subset_fastq_r1`, `subset_fastq_r2`, `subset_fastq_r3`: three fastq files containing 10000 
  reads aligned to chromosome 21 and 2000 reads aligned to other chromosomes or unaligned.  

# Limitations
- must be run single-threaded; use of additional cores with STAR causes output disordering,
  making index-based extraction as implemented here, impossible. 
