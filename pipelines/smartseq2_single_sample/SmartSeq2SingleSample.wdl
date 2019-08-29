import "HISAT2.wdl" as HISAT2
import "Picard.wdl" as Picard
import "RSEM.wdl" as RSEM
import "GroupMetricsOutputs.wdl" as GroupQCs
import "ZarrUtils.wdl" as ZarrUtils

workflow SmartSeq2SingleCell {
  meta {
    description: "Process SmartSeq2 scRNA-Seq data, include reads alignment, QC metrics collection, and gene expression quantitication"
  }
  # version of this pipeline
  String version = "smartseq2_v2.4.0"
  # load annotation
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat
  # load index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index
  # ref index name
  String hisat2_ref_name
  String hisat2_ref_trans_name
  # samples
  String stranded
  String sample_name
  String output_name
  File fastq1
  File fastq2

  # whether to convert the outputs to Zarr format, by default it's set to true
  Boolean output_zarr = true

  parameter_meta {
    genome_ref_fasta: "Genome reference in fasta format"
    rrna_intervals: "rRNA interval file required by Picard"
    gene_ref_flat: "Gene refflat file required by Picard"
    hisat2_ref_index: "HISAT2 reference index file in tarball"
    hisat2_ref_trans_index: "HISAT2 transcriptome index file in tarball"
    rsem_ref_index: "RSEM reference index file in tarball"
    hisat2_ref_name: "HISAT2 reference index name"
    hisat2_ref_trans_name: "HISAT2 transcriptome index file name"
    stranded: "Library strand information example values: FR RF NONE"
    sample_name: "Sample name or Cell ID"
    output_name: "Output name, can include path"
    fastq1: "R1 in paired end reads"
    fastq2: "R2 in paired end reads"
    output_zarr: "whether to run the taks that converts the outputs to Zarr format, by default it's true"
  }

  String quality_control_output_basename = output_name + "_qc"

  call HISAT2.HISAT2PairedEnd {
    input:
      hisat2_ref = hisat2_ref_index,
      fastq1 = fastq1,
      fastq2 = fastq2,
      ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_basename = quality_control_output_basename,
  }

  call Picard.CollectMultipleMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      genome_ref_fasta = genome_ref_fasta,
      output_basename = quality_control_output_basename,
  }

  call Picard.CollectRnaMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      ref_flat = gene_ref_flat,
      rrna_intervals = rrna_intervals,
      output_basename = quality_control_output_basename,
      stranded = stranded,
  }

  call Picard.CollectDuplicationMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      output_basename = quality_control_output_basename,
  }

  String data_output_basename = output_name + "_rsem"
  
  call HISAT2.HISAT2RSEM as HISAT2Transcriptome {
    input:
      hisat2_ref = hisat2_ref_trans_index,
      fastq1 = fastq1,
      fastq2 = fastq2,
      ref_name = hisat2_ref_trans_name,
      sample_name = sample_name,
      output_basename = data_output_basename,
  }

  call RSEM.RSEMExpression {
    input:
      trans_aligned_bam = HISAT2Transcriptome.output_bam,
      rsem_genome = rsem_ref_index,
      output_basename = data_output_basename,
      is_paired = true
  }

  call GroupQCs.GroupQCOutputs {
   input:
      picard_row_outputs = [CollectMultipleMetrics.alignment_summary_metrics,CollectMultipleMetrics.insert_size_metrics[0],CollectDuplicationMetrics.dedup_metrics,CollectRnaMetrics.rna_metrics,CollectMultipleMetrics.gc_bias_summary_metrics],
      picard_table_outputs = [CollectMultipleMetrics.base_call_dist_metrics,CollectMultipleMetrics.gc_bias_detail_metrics,CollectMultipleMetrics.pre_adapter_details_metrics,CollectMultipleMetrics.pre_adapter_summary_metrics,CollectMultipleMetrics.bait_bias_detail_metrics,CollectMultipleMetrics.error_summary_metrics],
      hisat2_stats = HISAT2PairedEnd.log_file,
      hisat2_trans_stats = HISAT2Transcriptome.log_file,
      rsem_stats = RSEMExpression.rsem_cnt,
      output_name = output_name
  }
  if (output_zarr) {
    call ZarrUtils.SmartSeq2ZarrConversion {
      input:
        rsem_gene_results = RSEMExpression.rsem_gene,
        smartseq_qc_files = GroupQCOutputs.group_files,
        sample_name=sample_name
    }
  }


  output {
    # version of this pipeline
    String pipeline_version = version
    # quality control outputs
    File aligned_bam = HISAT2PairedEnd.output_bam
    File bam_index = HISAT2PairedEnd.bam_index
    File insert_size_metrics = CollectMultipleMetrics.insert_size_metrics[0]
    File quality_distribution_metrics = CollectMultipleMetrics.quality_distribution_metrics
    File quality_by_cycle_metrics = CollectMultipleMetrics.quality_by_cycle_metrics
    File bait_bias_summary_metrics = CollectMultipleMetrics.bait_bias_summary_metrics
    File rna_metrics = CollectRnaMetrics.rna_metrics
    Array[File] group_results = GroupQCOutputs.group_files
    # data outputs
    File aligned_transcriptome_bam = HISAT2Transcriptome.output_bam
    File rsem_gene_results = RSEMExpression.rsem_gene
    File rsem_isoform_results = RSEMExpression.rsem_isoform

    # zarr
    Array[File]? zarr_output_files = SmartSeq2ZarrConversion.zarr_output_files
  }
}
