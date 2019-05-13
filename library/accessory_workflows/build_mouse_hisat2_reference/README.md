# Build HISAT2 Reference
Create a HISAT2 index to incorporate genome+snp+exon+splicing. HISAT2 provides script to build index with Ensembl or RefSeq index. We modify this script to incorporate with Gencode reference. 

# Input
## Genome Reference 
GRCh38 primary reference build, not include patch. 
## Gencode index
- gencode gtf version number. `27` is the current version.
- dbsnp version. current release is `150`
- output gencode reference bundle name. `hisat2_gencode_snp` 

## Ensembl index
- Ensembl gtf release number. `90` is the current release.
- dbsnp version. `150` is the current release.
- output Ensembl reference bundle name. `hisat2_ensembl_snp`
- 
