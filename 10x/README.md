## 10x Genomics Pipelines

This is code for running 10x data on the HCA Green Box. It perform the 
following tasks: 

1. Extracts barcode information from read 1 and index read 1
2. Aligns reads containing genomic information to a reference genome (using STAR)
3. Tags aligned bam reads with barcode information extracted in step 1
4. Identifies duplicate reads
5. Quantifies molecule counts in each gene and cell
6. Filters erroneous cell and molecular barcodes

# Download
Use git clone:
```
git clone git@github.com:HumanCellAtlas/skylab.git
cd skylab/10x/wdls/
```

# Requirements
## Dockers
- 10x clean docker: `marcusczi/cellranger_clean:latest`

## Inputs
- `read_paths`: path to a fastq folder containing all sequencing reads for the experiment
- `sample_def`: json file containing run information (see sample_def.json for example):  # TODO I omitted some here, should I have?
  - `sample_names`: name for experiment
  - `fastq_mode`: type of bcl demultiplexing to be done on the data
  - `read_paths`: location for fastq information #todo duplicate of above??
- `reads_per_file`: size of chunks (in reads) to create when scattering across processors
- `primers`: Illumina P5 and P7 sequencing adapters used
- `subsample_rate`: if passed and < 1.0, downsample data to this fraction. 
- `align`: The software to use to align reads
- `reference_path`: location of cellranger reference 
  (constructed with `skylab/10x/generate_reference_bundle.wdl`)

## Runtime Requirements
Memory: 50 GB  # todo is this right? 
Processors: 1+ (supports multiprocessing)
Disk Space: Input data size * 5 (normally approximately _ ) GB # todo marcus what's your estimate here? 

# Example Data
Example data was extracted from a public 10x human peripheral blood mononuclear cell 
experiment, where approximately 8000 cells were processed. Testing data was generated with 
`skylab/10x/generate_demo_data.wdl`, which extracts a small number of reads expected to align 
to chromosome 21 from the larger 10x experiment. A chromosome 21 reference bundle was created
with `skylab/10x/generate_reference_bundle.wdl`.    

# Output Description
## Run Summary HTML:
A high-level summary of the run results including analysis metrics 
## Gene-barcode matrices MEX and HDF5s
Count data is output in two formats: Sparse text files in matrix exchange format (.mtx) and HDF5 
database format. In both formats, the matrices are organized as barcodes (cells; rows) x genes 
(columns). In both formats, two files are provided: one containing all cells, and one filtered
to exclude those cell barcodes that failed analysis. 

## Barcoded BAMs
BAM formatted files including all analysis tags. see <a href=https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/bam>Barcoded BAMs</a>
for more information on specific tags. 

## Molecule info HDF5s
HDF databases containing per-molecule information for all molecules that contain a valid 
cell-barcode and valid UMI. The database schema is defined in more detail <a href=https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/output/molecule_info>here</a>
