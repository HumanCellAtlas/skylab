import "star.wdl" as star
import "picard.wdl" as picard
import "rsem.wdl" as rsem
import "featurecounts.wdl"  as featurecounts
import "htseq.wdl" as htseq

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
  Array[String] featureType = ['exon','gene','transcript']
	
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
  call picard.CollectMultipleMetrics {
    input:
      aligned_bam = Star.output_bam,
      ref_genome_fasta = ref_fasta,
      output_filename = output_prefix
  }
  call picard.CollectRnaMetrics {
    input:
      aligned_bam = Star.output_bam,
      ref_flat = ref_flat,
      rrna_interval = rrna_interval,
      output_filename = output_prefix,
      stranded = stranded
  }
  call picard.CollectDuplicationMetrics {
    input:
      aligned_bam = Star.output_bam,
      output_filename = output_prefix
  }
  call rsem.RsemExpression as Rsem {
    input:
      trans_aligned_bam = Star.output_bam_trans,
      rsem_genome = rsem_genome,
      rsem_out = output_prefix
  }
  call featurecounts.FeatureCountsUniqueMapping {
    input:
      aligned_bam = Star.output_bam,
      gtf = gtf, 
      fc_out = output_prefix
  }
  call featurecounts.FeatureCountsMultiMapping {
    input:
      aligned_bam = Star.output_bam,
      gtf = gtf,
      fc_out = output_prefix
  }
  scatter (ftype in featureType) {
    call htseq.htseq_count {
      input:
        aligned_bam = Star.output_bam,
        gtf = gtf,
        featuretype = ftype,
        output_filename = output_prefix
    }
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
    File gc_bias_detail_metrics =  CollectMultipleMetrics.gc_bias_detail_metrics
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
    File rna_coverage =CollectRnaMetrics.rna_coverage
    File dedup_metrics = CollectDuplicationMetrics.dedup_metrics
    File dedup_bamfile = CollectDuplicationMetrics.dedup_bamfile
    File unq_exons_counts = FeatureCountsUniqueMapping.exons
    File unq_genes_counts = FeatureCountsUniqueMapping.genes
    File unq_trans_counts = FeatureCountsUniqueMapping.trans
    File mult_exons_counts = FeatureCountsMultiMapping.exons
    File mult_genes_counts = FeatureCountsMultiMapping.genes
    File mult_trans_counts = FeatureCountsMultiMapping.trans
    File rsem_gene_results = Rsem.rsem_gene
    File rsem_isoform_results = Rsem.rsem_isoform
    File rsem_gene_counts = Rsem.rsem_gene_count
    File rsem_time_log = Rsem.rsem_time
    File rsem_cnt_log = Rsem.rsem_cnt
    File rsem_model_log = Rsem.rsem_model
    File rsem_theta_log = Rsem.rsem_theta
    Array[File] htseq_counts = htseq_count.counts
  }
}
