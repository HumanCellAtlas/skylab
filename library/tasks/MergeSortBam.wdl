task MergeSortBamFiles {
  Array[File] bam_inputs
  String sort_order

  Int disk_size = 500
  Int compression_level = 5

  command {
    set -e

    java -Dsamjdk.compression_level=${compression_level} -Xms7000m -jar /usr/gitc/picard.jar \
      MergeSamFiles \
      USE_THREADING=true \
      SORT_ORDER=${sort_order} \
      INPUT=${sep=' INPUT=' bam_inputs} \
      OUTPUT=merged.bam \
  }

  runtime {
    memory: "7.5 GB"
    cpu: 1
    disks: "local-disk ${disk_size} HDD"
    docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
  }
  output {
    File output_bam = "merged.bam"
  }
}
