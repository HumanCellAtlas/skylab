task CollectMultipleMetrics {
  File aligned_bam
  File genome_ref_fasta
  String output_basename

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"])
  Int machine_mem_mb = select_first([opt_memory_gb, 7]) * 1000
  # give the command 1 GB of overhead
  Int command_mem_mb = machine_mem_mb - 1000
  Int cpu = select_first([opt_cpu, 1])
  # use provided disk number or dynamically size on our own, with 10GB of additional disk
  Int disk = select_first([opt_disk, ceil(size(aligned_bam, "GB") + size(genome_ref_fasta, "GB") + 10)])
  Int preemptible = select_first([opt_preemptible, 5])

  meta {
    description: "JISHU HELP AGAINNNNNNNNNNNNN"
  }

  parameter_meta {
    aligned_bam: ""
    genome_ref_fasta: ""
    output_basename: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar CollectMultipleMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_basename}" \
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
      REFERENCE_SEQUENCE="${genome_ref_fasta}" \
      ASSUME_SORTED=true
  }
  
  runtime {
    docker: docker
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File alignment_summary_metrics = "${output_basename}.alignment_summary_metrics.txt"
    File base_call_dist_metrics = "${output_basename}.base_distribution_by_cycle_metrics.txt"
    File base_call_pdf = "${output_basename}.base_distribution_by_cycle.pdf"
    File gc_bias_detail_metrics = "${output_basename}.gc_bias.detail_metrics.txt"
    File gc_bias_dist_pdf = "${output_basename}.gc_bias.pdf"
    File gc_bias_summary_metrics = "${output_basename}.gc_bias.summary_metrics.txt"
    File insert_size_hist = "${output_basename}.insert_size_histogram.pdf"
    File insert_size_metrics = "${output_basename}.insert_size_metrics.txt"
    File quality_distribution_metrics = "${output_basename}.quality_distribution_metrics.txt"
    File quality_distribution_dist_pdf = "${output_basename}.quality_distribution.pdf"
    File quality_by_cycle_metrics = "${output_basename}.quality_by_cycle_metrics.txt"
    File quality_by_cycle_pdf = "${output_basename}.quality_by_cycle.pdf"
    File pre_adapter_details_metrics = "${output_basename}.pre_adapter_detail_metrics.txt"
    File pre_adapter_summary_metrics = "${output_basename}.pre_adapter_summary_metrics.txt"
    File bait_bias_detail_metrics = "${output_basename}.bait_bias_detail_metrics.txt"
    File bait_bias_summary_metrics = "${output_basename}.bait_bias_summary_metrics.txt"
    File error_summary_metrics = "${output_basename}.error_summary_metrics.txt"
  }
}

task CollectRnaMetrics {
  File aligned_bam
  File ref_flat
  File rrna_intervals
  String output_basename
  String stranded
  
  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  # give the command 500 MB of overhead
  Int command_mem_mb = machine_mem_mb - 500
  Int cpu = select_first([opt_cpu, 1])
  # use provided disk number or dynamically size on our own, with 10GB of additional disk
  Int disk = select_first([opt_disk, ceil(size(aligned_bam, "GB") + size(ref_flat, "GB") + size(rrna_intervals, "GB") + 10)])
  Int preemptible = select_first([opt_preemptible, 5])

  meta {
    description: "JISHU HELP AGAINNNNNNNNNNNNN"
  }

  parameter_meta {
    aligned_bam: ""
    ref_flat: ""
    rrna_intervals: ""
    output_basename: ""
    stranded: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }
  
  command {
    set -e
    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar CollectRnaSeqMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_basename}.rna_metrics.txt" \
      REF_FLAT="${ref_flat}" \
      RIBOSOMAL_INTERVALS="${rrna_intervals}" \
      STRAND_SPECIFICITY=${stranded} \
      CHART_OUTPUT="${output_basename}.rna.coverage.pdf"
    touch "${output_basename}.rna.coverage.pdf"
  }
  
  runtime {
    docker: docker
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File rna_metrics = "${output_basename}.rna_metrics.txt"
    File rna_coverage_pdf = "${output_basename}.rna.coverage.pdf"
  }
}

# Here we use "-XX:ParallelGCThreads=2" to run MarkDuplication on multiple threads 
task CollectDuplicationMetrics {
  File aligned_bam
  String output_basename

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"])
  Int machine_mem_mb = select_first([opt_memory_gb, 7]) * 1000
  # give the command 1 GB of overhead
  Int command_mem_mb = machine_mem_mb - 1000
  Int cpu = select_first([opt_cpu, 2])
  # use provided disk number or dynamically size on our own, with 10GB of additional disk
  Int disk = select_first([opt_disk, ceil(size(aligned_bam, "GB") + 10)])
  Int preemptible = select_first([opt_preemptible, 5])

  meta {
    description: "JISHU HELP AGAINNNNNNNNNNNNN"
  }

  parameter_meta {
    aligned_bam: ""
    output_basename: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }
  
  command {
    java -Xmx${command_mem_mb}m -XX:ParallelGCThreads=${cpu}  -jar /usr/picard/picard.jar  MarkDuplicates \
       VALIDATION_STRINGENCY=SILENT  \
       INPUT=${aligned_bam} \
       OUTPUT="${output_basename}.MarkDuplicated.bam" \
       ASSUME_SORTED=true \
       METRICS_FILE="${output_basename}.duplicate_metrics.txt" \
       REMOVE_DUPLICATES=false
  }
  
  runtime {
    docker: docker
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File dedup_metrics = "${output_basename}.duplicate_metrics.txt"
  }
}

