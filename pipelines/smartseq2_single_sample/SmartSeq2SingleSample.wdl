import "HISAT2.wdl" as HISAT2
import "Picard.wdl" as Picard
import "RSEM.wdl" as RSEM

workflow SmartSeq2SingleCell {
  meta {
    description: "Process SmartSeq2 scRNA-Seq data, include reads alignment, QC metrics collection, and gene expression quantitication"
  }

  # load annotation
  File gtf_file
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat

  # load index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index

  # ref index name
  File hisat2_ref_name
  File hisat2_ref_trans_name

  # samples
  String stranded
  String sample_name
  String output_name
  File fastq1
  File fastq2

  parameter_meta {
    gtf_file: ""
    genome_ref_fasta: ""
    rrna_intervals: ""
    gene_ref_flat: ""
    hisat2_ref_index: ""
    hisat2_ref_trans_index: ""
    rsem_ref_index: ""
    hisat2_ref_name: ""
    hisat2_ref_trans_name: ""
    stranded: ""
    sample_name: ""
    output_name: ""
    fastq1: ""
    fastq2: ""
  }

  String quality_control_output_prefix = output_name + "_qc"
  call HISAT2.HISAT2PairedEnd {
    input:
      hisat2_ref = hisat2_ref_index,
      fq1 = fastq1,
      fq2 = fastq2,
      ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_name = quality_control_output_prefix
  }

  call Picard.CollectMultipleMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      genome_ref_fasta = genome_ref_fasta,
      output_name = quality_control_output_prefix
  }

  call Picard.CollectRnaMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      ref_flat = gene_ref_flat,
      rrna_intervals = rrna_intervals,
      output_name = quality_control_output_prefix,
      stranded = stranded,
  }

  call Picard.CollectDuplicationMetrics {
    input:
      aligned_bam = HISAT2PairedEnd.output_bam,
      output_name = quality_control_output_prefix
  }

  String data_output_prefix = output_name + "_rsem"
  # TODO should `Trans` be Transcript??
  call HISAT2.HISAT2RSEM as HISAT2Trans {
    input:
      hisat2_ref = hisat2_ref_trans_index,
      fq1 = fastq1,
      fq2 = fastq2,
      ref_name = hisat2_ref_trans_name,
      sample_name = sample_name,
      output_name = data_output_prefix,
  }

  call RSEM.RSEMExpression {
    input:
      trans_aligned_bam = HISAT2Trans.output_bam,
      rsem_genome = rsem_ref_index,
      rsem_out = data_output_prefix,
  }

  output {
    # quality control outputs
    File aligned_bam = HISAT2PairedEnd.output_bam
    File met_file = HISAT2PairedEnd.met_file
    File log_file = HISAT2PairedEnd.log_file
    File alignment_summary_metrics = CollectMultipleMetrics.alignment_summary_metrics
    File base_call_dist_metrics = CollectMultipleMetrics.base_call_dist_metrics
    File base_call_pdf = CollectMultipleMetrics.base_call_pdf
    File gc_bias_detail_metrics = CollectMultipleMetrics.gc_bias_detail_metrics
    File gc_bias_dist_pdf = CollectMultipleMetrics.gc_bias_dist_pdf
    File gc_bias_summary_metrics = CollectMultipleMetrics.gc_bias_summary_metrics
    File insert_size_hist = CollectMultipleMetrics.insert_size_hist
    File insert_size_metrics = CollectMultipleMetrics.insert_size_metrics
    File quality_distribution_metrics = CollectMultipleMetrics.quality_distribution_metrics
    File quality_distribution_dist_pdf = CollectMultipleMetrics.quality_distribution_dist_pdf
    File quality_by_cycle_metrics = CollectMultipleMetrics.quality_by_cycle_metrics
    File quality_by_cycle_pdf = CollectMultipleMetrics.quality_by_cycle_pdf
    File pre_adapter_details_metrics = CollectMultipleMetrics.pre_adapter_details_metrics
    File bait_bias_detail_metrics = CollectMultipleMetrics.bait_bias_detail_metrics
    File bait_bias_summary_metrics = CollectMultipleMetrics.bait_bias_summary_metrics
    File error_summary_metrics = CollectMultipleMetrics.error_summary_metrics
    File rna_metrics = CollectRnaMetrics.rna_metrics
    File rna_coverage = CollectRnaMetrics.rna_coverage_pdf
    File dedup_metrics = CollectDuplicationMetrics.dedup_metrics

    # data outputs
    File aligned_trans_bam = HISAT2Trans.output_bam
    File hisat2tran_met_file = HISAT2Trans.met_file
    File hisat2tran_log_file = HISAT2Trans.log_file
    File rsem_gene_results = RSEMExpression.rsem_gene
    File rsem_isoform_results = RSEMExpression.rsem_isoform
    File rsem_time_log = RSEMExpression.rsem_time
    File rsem_cnt_log = RSEMExpression.rsem_cnt
    File rsem_model_log = RSEMExpression.rsem_model
    File rsem_theta_log = RSEMExpression.rsem_theta
  }
}
