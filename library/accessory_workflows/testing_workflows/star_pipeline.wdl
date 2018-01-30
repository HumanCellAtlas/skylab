import "star.wdl" as star
import "picard.wdl" as picard
import "rsem.wdl" as rsem
import "featurecounts.wdl"  as featurecounts

workflow RunStarPipeline {
  File fastq_read1
  File fastq_read2
  File gtf
  File stranded
  File ref_fasta
  File rrna_interval
  File ref_flat 
  File star_genome
  File rsem_genome
  String sample_name
  String output_prefix
	## estimate size
  Float star_ref_size = size(star_genome, "GB")
  Float fastq_size = size(fastq_read1, "GB") +size(fastq_read2, "GB")
  Float rsem_ref_size = size(rsem_genome, "GB")
  Float reference_bundle_size = size(ref_fasta, "GB") + size(ref_flat, "GB") + size(rrna_interval, "GB") + size(gtf, "GB")
  Float bam_disk_multiplier = 10.0
  Int? increase_disk_size
  Int additional_disk = select_first([increase_disk_size, 10])

  call star.StarPE as Star {
    input:
      input_fastq_read1 = fastq_read1,
      input_fastq_read2 = fastq_read2,
      gtf = gtf,
      star_genome = star_genome,
      sample_tag = output_prefix,
      pu_tag = sample_name,
      lib_tag = sample_name,
      id_tag = sample_name
  }
  Float bam_size = size(Star.output_bam, "GB") 
  call picard.CollectMultipleMetrics {
    input:
      aligned_bam = Star.output_bam,
      ref_genome_fasta = ref_fasta,
      output_filename = output_prefix,
      disk_size = ceil(bam_size + reference_bundle_size + additional_disk)

  }
  call picard.CollectRnaMetrics {
    input:
      aligned_bam = Star.output_bam,
      ref_flat = ref_flat,
      rrna_interval = rrna_interval,
      output_filename = output_prefix,
      stranded = stranded,
      disk_size = ceil(bam_size + reference_bundle_size + additional_disk)
  }
  call picard.CollectDuplicationMetrics {
    input:
      aligned_bam = Star.output_bam,
      output_filename = output_prefix,
      disk_size = ceil((bam_disk_multiplier * bam_size) + additional_disk)
  }
  call rsem.RsemExpression as Rsem {
    input:
      trans_aligned_bam = Star.output_bam_trans,
      rsem_genome = rsem_genome,
      rsem_out = output_prefix,
      disk_size = ceil(fastq_size * bam_disk_multiplier+rsem_ref_size+additional_disk * 2.0)
  }
  call featurecounts.FeatureCountsUniqueMapping {
    input:
      aligned_bam = Star.output_bam,
      gtf = gtf, 
      fc_out = output_prefix,
      disk_size = ceil((bam_disk_multiplier * bam_size) + additional_disk * 2.0)
  }
  call featurecounts.FeatureCountsMultiMapping {
    input:
      aligned_bam = Star.output_bam,
      gtf = gtf,
      fc_out = output_prefix,
      disk_size = ceil((bam_disk_multiplier * bam_size) + additional_disk * 2.0)
  }
  output {
    File aligned_bam = Star.output_bam
    File aligned_trans_bam = Star.output_bam_trans
    File junction_table = Star.junction_table
    File final_log = Star.final_log
    File process_log = Star.process_log
    File logfile = Star.logfile
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
    File unq_exons_counts = FeatureCountsUniqueMapping.exons
    File unq_genes_counts = FeatureCountsUniqueMapping.genes
    File unq_trans_counts = FeatureCountsUniqueMapping.trans
    File mult_exons_counts = FeatureCountsMultiMapping.exons
    File mult_genes_counts = FeatureCountsMultiMapping.genes
    File mult_trans_counts = FeatureCountsMultiMapping.trans
    File rsem_gene_results = Rsem.rsem_gene
    File rsem_isoform_results = Rsem.rsem_isoform
    File rsem_time_log = Rsem.rsem_time
    File rsem_cnt_log = Rsem.rsem_cnt
    File rsem_model_log = Rsem.rsem_model
    File rsem_theta_log = Rsem.rsem_theta
  }
}
