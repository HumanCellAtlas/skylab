task CollectMultipleMetrics {
  File aligned_bam
  File ref_genome_fasta
  String output_filename
  Int disk_size 
  command {
    java -Xmx6g -jar /usr/picard/picard.jar CollectMultipleMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_filename}" \
      FILE_EXTENSION=".txt" \
      PROGRAM=null \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectGcBiasMetrics \
      PROGRAM=CollectBaseDistributionByCycle \
      PROGRAM=QualityScoreDistribution \
      PROGRAM=MeanQualityByCycle \
      PROGRAM=CollectSequencingArtifactMetrics \
      PROGRAM=CollectQualityYieldMetrics \
      REFERENCE_SEQUENCE="${ref_genome_fasta}" \
      ASSUME_SORTED=true
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    memory:"7.5 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }
  output {
    File alignment_summary_metrics = "${output_filename}.alignment_summary_metrics.txt"
    File base_call_dist_metrics = "${output_filename}.base_distribution_by_cycle_metrics.txt"
    File base_call_pdf = "${output_filename}.base_distribution_by_cycle.pdf"
    File gc_bias_detail_metrics = "${output_filename}.gc_bias.detail_metrics.txt"
    File gc_bias_dist_pdf = "${output_filename}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${output_filename}.gc_bias.summary_metrics.txt"
    File insert_size_hist = "${output_filename}.insert_size_histogram.pdf"
    File insert_size_metrics = "${output_filename}.insert_size_metrics.txt"
    File quality_distribution_metrics = "${output_filename}.quality_distribution_metrics.txt"
    File quality_distribution_dist_pdf = "${output_filename}.quality_distribution.pdf"
    File quality_by_cycle_metrics = "${output_filename}.quality_by_cycle_metrics.txt"
    File quality_by_cycle_pdf = "${output_filename}.quality_by_cycle.pdf"
    File pre_adapter_details_metrics = "${output_filename}.pre_adapter_detail_metrics.txt"
    File pre_adapter_summary_metrics = "${output_filename}.pre_adapter_summary_metrics.txt"
    File bait_bias_detail_metrics = "${output_filename}.bait_bias_detail_metrics.txt"
    File bait_bias_summary_metrics = "${output_filename}.bait_bias_summary_metrics.txt"
    File error_summary_metrics = "${output_filename}.error_summary_metrics.txt"
  }
}

task CollectRnaMetrics {
  File aligned_bam
  File ref_flat
  File rrna_interval
  String output_filename
  String stranded
  Int disk_size
  command{
    set -e
    java -Xmx3g -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_filename}.rna_metrics.txt" \
      REF_FLAT="${ref_flat}" \
      RIBOSOMAL_INTERVALS="${rrna_interval}" \
      STRAND_SPECIFICITY=${stranded} \
      CHART_OUTPUT="${output_filename}.rna.coverage.pdf"
    touch "${output_filename}.rna.coverage.pdf"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    memory:"3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }
  output {
    File rna_metrics = "${output_filename}.rna_metrics.txt"
    File rna_coverage_pdf = "${output_filename}.rna.coverage.pdf"
  }
}
## Here are use  -XX:ParallelGCThreads=2 to run MarkDuplication on mutlple
## thread. 
task CollectDuplicationMetrics {
  File aligned_bam
  String output_filename
  Int disk_size
  command {
    java -Xmx6g -XX:ParallelGCThreads=2  -jar /usr/picard/picard.jar  MarkDuplicates \
       VALIDATION_STRINGENCY=SILENT  \
       INPUT=${aligned_bam} \
       OUTPUT="${output_filename}.MarkDuplicated.bam" \
       ASSUME_SORTED=true \
       METRICS_FILE="${output_filename}.duplicate_metrics.txt" \
       REMOVE_DUPLICATES=false
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    memory: "7.5 GB"
    cpu: "2"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 5
  }
  output {
    File dedup_metrics = "${output_filename}.duplicate_metrics.txt"
  }
}

