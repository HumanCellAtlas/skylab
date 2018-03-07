import "HISAT2.wdl" as HISAT2
import "Picard.wdl" as Picard

workflow RunHisat2Pipeline {
  meta {
    description: "JISHU HELP AGAINNNNNNNNNNNNN"
  }

  File fastq_read1
  File fastq_read2
  File gtf
  File stranded
  File ref_fasta
  File rrna_interval
  File ref_flat
  File hisat2_ref
  String output_prefix
  String hisat2_ref_name
  String sample_name

  parameter_meta {
    fastq_read1: ""
    fastq_read2: ""
    gtf: ""
    stranded: ""
    ref_fasta: ""
    rrna_interval: ""
    ref_flat: ""
    hisat2_ref: ""
    output_prefix: ""
    hisat2_ref_name: ""
    sample_name: ""
  }
 
  call HISAT2.HISAT2PairedEnd {
    input:
      hisat2_ref = hisat2_ref,
      fq1 = fastq_read1,
      fq2 = fastq_read2,
      ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_name = output_prefix
  }
  
  call Picard.CollectMultipleMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      ref_genome_fasta = ref_fasta,
      output_filename = output_prefix,
      disk_size = ceil(bam_size + reference_bundle_size + additional_disk)
  }

  call Picard.CollectRnaMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      ref_flat = ref_flat,
      rrna_interval = rrna_interval,
      output_filename = output_prefix,
      stranded = stranded,
      disk_size = ceil(bam_size + reference_bundle_size + additional_disk)
  }

  call Picard.CollectDuplicationMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      output_filename = output_prefix,
      disk_size = ceil((bam_disk_multiplier * bam_size) + additional_disk)
  }

  output {
    File aligned_bam = HISAT2PairedEnd.output_bam
    File metfile = HISAT2PairedEnd.metfile
    File logfile = HISAT2PairedEnd.logfile
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
  }
}
