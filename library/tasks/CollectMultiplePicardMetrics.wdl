
task CollectMultipleMetrics {
  File aligned_bam
  File ref_genome_fasta
  String output_filename
  Int estimated_required_disks = ceil((size(aligned_bam, "G") + size(ref_genome_fasta, "G")) * 1.2)

  command {
    java -Xmx3g -jar /usr/picard/picard.jar CollectMultipleMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT="${aligned_bam}" \
      OUTPUT="${output_filename}" \
      PROGRAM=CollectAlignmentSummaryMetrics \
      PROGRAM=CollectInsertSizeMetrics \
      PROGRAM=CollectGcBiasMetrics \
      REFERENCE_SEQUENCE="${ref_genome_fasta}"

      tar -czvf "${output_filename}.tar.gz" "${output_filename}"*
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    memory:"3.75 GB"
    disks: "local-disk ${estimated_required_disks} HDD"
  }
  output {
    File alignment_metrics = "${output_filename}.tar.gz"
  }
}
