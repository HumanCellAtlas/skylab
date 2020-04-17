version 1.0

import "SmartSeq2SingleSample.wdl" as single_cell_run
import "SmartSeq2PlateAggregation.wdl" as ss2_plate_aggregation
import "ZarrUtils.wdl" as ZarrUtils
       
workflow MultiSampleSmartSeq2 {
  meta {
    description: "The MultiSampleSmartSeq2 pipeline runs multiple SS2 samples in a single pipeline invocation"
  }

  input {
      # Version of this pipeline
      String version = "MultiSampleSmartSeq2_v1.0.0"

      # Gene Annotation
      File genome_ref_fasta
      File rrna_intervals
      File gene_ref_flat

      # Reference index information
      File hisat2_ref_name
      File hisat2_ref_trans_name
      File hisat2_ref_index
      File hisat2_ref_trans_index
      File rsem_ref_index

      # Sample information
      String stranded
      String file_prefix
      Array[File] read1
      Array[File]? read2
      String batch_id
      Boolean paired_end

      Boolean output_loom
  }

  # Parameter metadata information
  parameter_meta {
    genome_ref_fasta: "Genome reference in fasta format"
    rrna_intervals: "rRNA interval file required by Picard"
    gene_ref_flat: "Gene refflat file required by Picard"
    hisat2_ref_name: "HISAT2 reference index name"
    hisat2_ref_trans_name: "HISAT2 transcriptome index file name"
    hisat2_ref_index: "HISAT2 reference index file in tarball"
    hisat2_ref_trans_index: "HISAT2 transcriptome index file in tarball"
    rsem_ref_index: "RSEM reference index file in tarball"
    stranded: "Library strand information example values: FR RF NONE"
    fastq1: "Array of read1 fastq files"
    fastq2: "Optional array of read2 fastq files. "
    batch_id: " Identifier for the batch"
    paired_end: "Is the sample paired end or not"
  }

  ### Execution starts here ###
  if (paired_end) {
      scatter(idx in range(length(read1))) {
        call single_cell_run.SmartSeq2SingleCell as sc_pe {
          input:
            fastq1 = read1[idx],
            fastq2 = read2[idx],
            stranded = stranded,
            genome_ref_fasta = genome_ref_fasta,
            rrna_intervals = rrna_intervals,
            gene_ref_flat = gene_ref_flat,
            hisat2_ref_index = hisat2_ref_index,
            hisat2_ref_name = hisat2_ref_name,
            hisat2_ref_trans_index = hisat2_ref_trans_index,
            hisat2_ref_trans_name = hisat2_ref_trans_name,
            rsem_ref_index = rsem_ref_index,
            sample_name = basename(read1[idx], ".fastq.gz"),
            output_name = basename(read1[idx], ".fastq.gz"),
            paired_end = paired_end,
            output_zarr = true
        }
      }
  }
  if (!paired_end) {
        scatter(idx in range(length(read1))) {
          call single_cell_run.SmartSeq2SingleCell as sc_se {
            input:
              fastq1 = read1[1],
              stranded = stranded,
              genome_ref_fasta = genome_ref_fasta,
              rrna_intervals = rrna_intervals,
              gene_ref_flat = gene_ref_flat,
              hisat2_ref_index = hisat2_ref_index,
              hisat2_ref_name = hisat2_ref_name,
              hisat2_ref_trans_index = hisat2_ref_trans_index,
              hisat2_ref_trans_name = hisat2_ref_trans_name,
              rsem_ref_index = rsem_ref_index,
              sample_name = basename(read1[idx], ".fastq.gz"),
              output_name = basename(read1[idx], ".fastq.gz"),
              paired_end = paired_end,
              output_zarr = true
          }
        }
  }

  Array[Array[File]?] zarr_output_files = select_first([sc_pe.zarr_output_files, sc_se.zarr_output_files])
  Array[File] bam_files_intermediate = select_first([sc_pe.aligned_bam, sc_se.aligned_bam])
  Array[File] bam_index_files_intermediate = select_first([sc_pe.bam_index, sc_se.bam_index])

  ### Aggregate the Zarr Files Directly ###
  call ss2_plate_aggregation.AggregateSmartSeq2Zarr as AggregateZarr {
    input:
      zarr_input = zarr_output_files,
      output_file_name = batch_id
  }

   if (output_loom) {
    call ZarrUtils.SmartSeq2PlateToLoom as ZarrToLoom {
       input:
         batch_id = batch_id,
         zarr_files = AggregateZarr.zarr_output_files
    }
  }

  ### Pipeline output ###
  output {
    # Bam files and their indexes
    Array[File] bam_files = bam_files_intermediate
    Array[File] bam_index_files = bam_index_files_intermediate
    Array[File] zarrout = AggregateZarr.zarr_output_files
    File? loom_output = ZarrToLoom.loom_output
  }
}
