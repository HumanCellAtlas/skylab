# 10x/generate_reference_bundle 0.0.1
This pipeline performs the following tasks:
- Task 1 (using 10xGenomics/cellranger), calls cellranger mkref to generate a reference bundle for
use with cellranger count. It countains a STAR index and several additional 10x specific files. 

## Download
Use git clone: 

```
git clone git@github.com:HumanCellAtlas/skylab.git
cd skylab/10x/generate_reference_bundle
```

# Requirements
## Dockers
- cellranger docker (latest): `marcusczi/cellranger_clean:latest`, a docker containing public 10x 
  code. 

## File Inputs
- `fasta`: genome fasta file
- `gtf`: gtf file for matching organism
- `reference_name`: a name tag for the reference being generated

## Runtime requirements
- Memory: 60GB
- Processors: 1+ (multiprocessing supported)
- Disk Space: ~ 50GB
- Expected time: ~ 5 minutes

# Example Input Data
Input data was generated using wdls from `skylab` which subset  # todo which ones?
only chromosome 21 from the human genome and annotation references downloaded from 
<a href=https://www.gencodegenes.org/releases/current.html>GENCODE</a>

# Output Description
This WDL produces a tarred reference bundle for use with 10x pipelines. 
- `${reference_name}.tar` tarred 10x reference bundle
