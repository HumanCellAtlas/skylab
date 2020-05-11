version 1.0

import "https://raw.githubusercontent.com/HumanCellAtlas/skylab/jw_MultiSampleSmartSeq2_Terra_arrays/pipelines/smartseq2_single_sample/SmartSeq2SingleSample.wdl" as single_cell_run
import "https://raw.githubusercontent.com/HumanCellAtlas/skylab/jw_MultiSampleSmartSeq2_Terra_arrays/library/tasks/SmartSeq2PlateAggregation.wdl" as ss2_plate_aggregation
import "https://raw.githubusercontent.com/HumanCellAtlas/skylab/jw_MultiSampleSmartSeq2_Terra_arrays/library/tasks/ZarrUtils.wdl" as ZarrUtils
       
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
      Array[String] sample_names
      Array[String] fastq1_input_files
      Array[String] fastq2_input_files = []
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
    sample_names: "Array of sample names"
    fastq1_input_files: "Array of fastq1 files; order must match the order in sample_names."
    fastq2_input_files: "Array of fastq2 files for paired end runs; order must match fastq1_input_files and sample_names."
    batch_id: " Identifier for the batch"
    paired_end: "Is the sample paired end or not"
  }

  # Check that all input arrays are the same length
  call checkInputArrays as checkArrays{
      input:
         sample_names = sample_names,
         fastq1_input_files = fastq1_input_files,
         fastq2_input_files = fastq2_input_files
  }

  ### Execution starts here ###
  if (paired_end) {
      scatter(idx in range(length(sample_names))) {
        call single_cell_run.SmartSeq2SingleCell as sc_pe {
          input:
            fastq1 = fastq1_input_files[idx],
            fastq2 = fastq2_input_files[idx],
            stranded = stranded,
            genome_ref_fasta = genome_ref_fasta,
            rrna_intervals = rrna_intervals,
            gene_ref_flat = gene_ref_flat,
            hisat2_ref_index = hisat2_ref_index,
            hisat2_ref_name = hisat2_ref_name,
            hisat2_ref_trans_index = hisat2_ref_trans_index,
            hisat2_ref_trans_name = hisat2_ref_trans_name,
            rsem_ref_index = rsem_ref_index,
            sample_name = sample_names[idx],
            output_name = sample_names[idx],
            paired_end = paired_end,
            output_zarr = true
        }
      }
  }
  if (!paired_end) {
        scatter(idx in range(length(sample_names))) {
          call single_cell_run.SmartSeq2SingleCell as sc_se {
            input:
              fastq1 = fastq1_input_files[idx],
              stranded = stranded,
              genome_ref_fasta = genome_ref_fasta,
              rrna_intervals = rrna_intervals,
              gene_ref_flat = gene_ref_flat,
              hisat2_ref_index = hisat2_ref_index,
              hisat2_ref_name = hisat2_ref_name,
              hisat2_ref_trans_index = hisat2_ref_trans_index,
              hisat2_ref_trans_name = hisat2_ref_trans_name,
              rsem_ref_index = rsem_ref_index,
              sample_name = sample_names[idx],
              output_name = sample_names[idx],
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

task checkInputArrays {
  input {
    Array[String] sample_names
    Array[String] fastq1_input_files
    Array[String] fastq2_input_files
  }
  Int len_sample_names = length(sample_names)
  Int len_fastq1_input_files = length(fastq1_input_files)
  Int len_fastq2_input_files = length(fastq2_input_files)

  meta {
    description: "checks input arrays to ensure that all arrays are the same length"
  }

  command {
    set -e

    if [[ ~{len_sample_names} !=  ~{len_fastq1_input_files} ]]
      then
      echo "ERROR: Different number of arguments for sample names and fastq1 files"
      exit 1;
    fi

    if [[ ~{len_fastq2_input_files} != 0 && ~{len_fastq2_input_files} != ~{len_sample_names} ]]
      then
      echo "ERROR: Different number of arguments for sample names and fastq1 files"
      exit 1;
    fi
    exit 0;
  }

  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1 GiB"
    disks: "local-disk 1 HDD"
  }

}