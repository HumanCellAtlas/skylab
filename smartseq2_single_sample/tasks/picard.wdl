task CollectMultipleMetrics {
  File aligned_bam
  File ref_genome_fasta
  String output_filename
  
  command {
    java -Xmx4g -jar /usr/picard/picard.jar CollectMultipleMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT=${aligned_bam} \
      OUTPUT="${output_filename}" \
      FILE_EXTENSION = ".txt" \
      PROGRAM = CollectAlignmentSummaryMetrics \
      PROGRAM = CollectInsertSizeMetrics \
      PROGRAM = CollectGcBiasMetrics \
      PROGRAM = RnaSeqMetrics \
      REFERENCE_SEQUENCE=${ref_genome_fasta} \
      ASSUME_SORTED=true

      tar -cvf "${output_filename}.tar" ${output_filename}.*.txt
  }
  runtime {
    docker:"humancellatlas/picard:2.10.10"
    memory:"3.75 GB"
    disks: "local-disk 10 HDD"
  }
  output {
    File alignment_metrics = "${output_filename}.tar"
  }
}


task CollectDuplicationMetrics {
  File aligned_bam
  String output_filename

  command {
    java -Xmx4g -jar picard.jar  MarkDuplicates \
       VALIDATION_STRINGENCY=SILENT  \
       INPUT=${aligned_bam} \
       OUTPUT="${output_filename}.MarkDuplicated.bam" \
       ASSUME_SORTED=true \
       METRICS_FILE="${output_filename}.duplicate_metrics.txt" \
       REMOVE_DUPLICATES=false
  }
  runtime {
    docker: "humancellatlas/picard:2.10.10"
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
  }
  output {
    File dedup_metrics = "${output_filename}.duplicate_metrics.txt"
    File dedup_bamfile = "${output_filename}.MarkDuplicated.bam"
  }
}

