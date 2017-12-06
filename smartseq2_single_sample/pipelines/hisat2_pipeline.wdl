import "hisat2.wdl" as hisat2
import "picard.wdl" as picard
import "featurecounts.wdl"  as featurecounts
import "htseq.wdl" as htseq

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
  Array[String] featuretype = ['exon','gene','transcript']
  
  call hisat2.HISAT2PE as Hisat2 {
    input:
      hisat2_ref = hisat2_ref,
      fq1 = fastq_read1,
      fq2 = fastq_read2,
      ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_name = output_prefix
  }
  call picard.CollectMultipleMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      ref_genome_fasta = ref_fasta,
      output_filename = output_prefix
  }
  call picard.CollectRnaMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      ref_flat = ref_flat,
      rrna_interval = rrna_interval,
      output_filename = output_prefix,
      stranded = stranded
  }
  call picard.CollectDuplicationMetrics {
    input:
      aligned_bam = Hisat2.output_bam,
      output_filename = output_prefix
  }
  call featurecounts.FeatureCountsUniqueMapping {
    input:
      aligned_bam = Hisat2.output_bam,
      gtf = gtf, 
      fc_out = output_prefix
  }
  call featurecounts.FeatureCountsMultiMapping {
    input:
      aligned_bam = Hisat2.output_bam,
      gtf = gtf,
      fc_out = output_prefix
  }
  scatter (ftype in featuretype) {
    call htseq.htseq_count {
      input:
        aligned_bam = Hisat2.output_bam,
        gtf = gtf,
        featuretype = ftype,
        output_filename = output_prefix
    }
  }
  output {
    File aligned_bam = Hisat2.output_bam
    File metfile = Hisat2.metfile
    File logfile = Hisat2.logfile
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
    Array[File] htseq_counts = htseq_count.counts
  }
}
