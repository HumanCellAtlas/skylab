
task FastqToUBam {
  File gz_fastq_file  # input fastq file
  String sample_name  # name of sample matching this file, inserted into read group header

  # estimate that bam is approximately equal in size to fastq, add 20% buffer
  Int estimated_disk_required = ceil(size(gz_fastq_file, "G") * 2.2)

  command {
    java -Xmx2g -jar picard.jar FastqToSam \
      FASTQ="${gz_fastq_file}" \
      SORT_ORDER=unsorted \
      OUTPUT=bamfile.bam \
      SAMPLE_NAME="${sample_name}"
  }
  
  runtime {
    docker: "humancellatlas/picard:2.10.10"
    cpu: 1
    memory: "2.5 GB"
    disks: "local-disk ${estimated_disk_required} HDD"
  }
  
  output {
    File bam = "bamfile.bam"
  }
}
