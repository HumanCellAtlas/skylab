import "hisat2.wdl" as hisat2
import "picard.wdl" as picard
## QC pipeline
## HISAT2 as aligner to align reads to genome reference
## output: genome reference aligned bam file 
## Picard will produce a set of QC metricss
## output: a set of metrics files. 
workflow RunHisat2Pipeline {
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
  ## variables to estimate disk size
  Float hisat2_ref_size = size(hisat2_ref, "GB")
  Float fastq_size = size(fastq_read1, "GB") + size(fastq_read2, "GB")
  Float reference_bundle_size = size(ref_fasta, "GB") + size(ref_flat, "GB") + size(rrna_interval, "GB") + size(gtf, "GB")
  Float bam_disk_multiplier = 10.0
  Int? increase_disk_size
  Int additional_disk = select_first([increase_disk_size, 10])
 
  call hisat2.HISAT2PE as Hisat2 {
    input:
      hisat2_ref = hisat2_ref,
      fq1 = fastq_read1,
      fq2 = fastq_read2,
      ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_name = output_prefix,
      disk_size = ceil(hisat2_ref_size + fastq_size * bam_disk_multiplier + additional_disk * 5.0)
  }
  
  Float bam_size = size(Hisat2.output_bam, "GB")
  
  call picard.CollectMultipleMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      ref_genome_fasta = ref_fasta,
      output_filename = output_prefix,
      disk_size = ceil(bam_size + reference_bundle_size + additional_disk)
  }
  call picard.CollectRnaMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      ref_flat = ref_flat,
      rrna_interval = rrna_interval,
      output_filename = output_prefix,
      stranded = stranded,
      disk_size = ceil(bam_size + reference_bundle_size + additional_disk)
  }
  call picard.CollectDuplicationMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      output_filename = output_prefix,
      disk_size = ceil((bam_disk_multiplier * bam_size) + additional_disk)
  }

  output {
    File aligned_bam = Hisat2.output_bam
    File metfile = Hisat2.metfile
    File logfile = Hisat2.logfile
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
