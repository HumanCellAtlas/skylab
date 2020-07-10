# What's in the Smart-seq2 Multi Sample Pipeline Loom File?

The Loom file is an HDF5 file generated using [Loompy v.3.0.6](http://loompy.org/). It contains global attributes containing the plateid ([Table 1](#table-1-global-attributes)). The main matrices contain RSEM expected_counts and RSEM TPMs. The file additionally contains multiple metrics for both individual cells (the columns of the matrix; [Table 2](#table-2-column-attributes-cell-metrics)) and individual genes (the rows of the matrix; [Table 3](#table-3-row-attributes-gene-metrics)). The tables below document these metrics, list which tools generate them, and define them. This Loom file is the default matrix output of the Smart-seq2 Multi Sample pipeline.  

**Note**: Loom files generated by the Smart-seq2 Multi Sample are different from the final Loom file distributed on the [Human Cell Atlas Data Portal](https://data.humancellatlas.org/explore/projects), which removes some of the metadata detailed in this document and contains additional metadata relating to each individual project. 

## Table 1. Global Attributes
The global attributes in the Loom apply to the whole file, not any specific part. There are three global attributes for the Multip Sample Smart-seq2 Loom. 

| Attribute | Details |
| :-- | :-- |
| LOOM_SPEC_VERSION | String with the loom file spec version |
| CreationDate | Date Loom file was generated |
| sample_id | The plateid listed in the pipeline configuration file; by default it is set to the batch_id. |
 

## Table 2. Column Attributes (Cell Metrics) 

| Cell Metrics | Program |Details |
|:---|:---:|:---| 
| `ACCUMULATION_LEVEL` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `ALIGNED_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `AT_DROPOUT` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `Aligned 0 time` | [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml) | Number and percent reads aligned 0 times |
| `Aligned 1 time` | [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml) | Number and percent reads aligned 1 time |
| `Aligned >1 times` | [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml) | Number and percent reads aligned more than 1 time |
| `BAD_CYCLES.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `CODING_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `CORRECT_STRAND_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `ESTIMATED_LIBRARY_SIZE` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `GC_DROPOUT` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `GC_NC_0_19` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `GC_NC_20_39` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `GC_NC_40_59` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `GC_NC_60_79` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `GC_NC_80_100` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `IGNORED_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `INCORRECT_STRAND_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `INTERGENIC_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `INTRONIC_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `MEAN_READ_LENGTH.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `MEDIAN_3PRIME_BIAS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `MEDIAN_5PRIME_BIAS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `MEDIAN_5PRIME_TO_3PRIME_BIAS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `MEDIAN_CV_COVERAGE` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `NUM_R1_TRANSCRIPT_STRAND_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `NUM_R2_TRANSCRIPT_STRAND_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `NUM_UNEXPLAINED_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `Overall alignment rate` | [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml) | Overall percent of reads that aligned |
| `PCT_ADAPTER.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_CHIMERAS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_CODING_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_CORRECT_STRAND_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_INTERGENIC_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_INTRONIC_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_MRNA_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_PF_READS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_PF_READS_ALIGNED.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_PF_READS_IMPROPER_PAIRS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_R1_TRANSCRIPT_STRAND_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_R2_TRANSCRIPT_STRAND_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_READS_ALIGNED_IN_PAIRS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_RIBOSOMAL_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_USABLE_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PCT_UTR_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PERCENT_DUPLICATION` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_ALIGNED_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_ALIGNED_BASES.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_HQ_ALIGNED_BASES.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_HQ_ALIGNED_Q20_BASES.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_HQ_ALIGNED_READS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_HQ_ERROR_RATE.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_HQ_MEDIAN_MISMATCHES.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_INDEL_RATE.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_MISMATCH_RATE.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_NOISE_READS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_READS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_READS_ALIGNED.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `PF_READS_IMPROPER_PAIRS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `READS_ALIGNED_IN_PAIRS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `READS_USED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `READ_PAIRS_EXAMINED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `READ_PAIR_DUPLICATES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `READ_PAIR_OPTICAL_DUPLICATES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `RIBOSOMAL_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `SECONDARY_OR_SUPPLEMENTARY_RDS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `STRAND_BALANCE.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `TOTAL_CLUSTERS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `TOTAL_READS.UNPAIRED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `Total reads` | [HISAT2](https://ccb.jhu.edu/software/hisat2/manual.shtml) | Total number of aligned reads |
| `UNMAPPED_READS` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `UNPAIRED_READS_EXAMINED` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `UNPAIRED_READ_DUPLICATES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `UTR_BASES` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `WINDOW_SIZE` | [Picard](https://broadinstitute.github.io/picard/picard-metric-definitions.html) |
| `alignable reads` |
| `cell_names` |
| `filtered reads` |
| `multiple mapped` |
| `strand` |
| `total alignments` |
| `total reads` |
| `unalignable reads` |
| `uncertain reads` |
| `unique aligned` |


## Table 3. Row Attributes (Gene Metrics)

| Gene Metrics                  | Program            |Details                 | 
|-------------------------------|--------------------|------------------------|
|`ensembl_ids` | [GENCODE GTF](https://www.gencodegenes.org/) | The gene_id listed in the GENCODE GTF |
|`gene_names` | [GENCODE GTF](https://www.gencodegenes.org/) | The unique gene_name provided in the GENCODE GTF |
