
task FastqToUBam {
  File fastq_file  # input fastq file
  String sample_id  # name of sample matching this file, inserted into read group header

  # a suffix to add to the fastq file; useful with mangled file IDs, since picard requires that
  # the file end in .gz or it will not detect the gzipping.
  String fastq_suffix = ""

  # estimate that bam is approximately equal in size to fastq, add 20% buffer
  Int estimated_disk_required = ceil(size(fastq_file, "G") * 2.2)

  command {
    mv "${fastq_file}" "${fastq_file}""${fastq_suffix}"  # add suffix; does nothing if not provided
    java -Xmx2g -jar /usr/picard/picard.jar FastqToSam \
      FASTQ="${fastq_file}""${fastq_suffix}" \
      SORT_ORDER=unsorted \
      OUTPUT=bamfile.bam \
      SAMPLE_NAME="${sample_id}"
  }
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    cpu: 1
    memory: "2.5 GB"
    disks: "local-disk ${estimated_disk_required} HDD"
  }
  
  output {
    File bam_output = "bamfile.bam"
  }
}
