
task MergeSortBamFiles {
  Array[Array[File]] bam_inputs
  Int disk_size = 500
  Int compression_level = 5

  command {
    input_line=$(python -c "print(' INPUT='.join('${sep=' ' bam_inputs}'.replace('[', '').replace(']', '').replace(',', '').split()))")

    echo $input_line
    java -Dsamjdk.compression_level=${compression_level} -Xms7000m -jar /usr/gitc/picard.jar \
      MergeSamFiles \
      USE_THREADING=true \
      SORT_ORDER=coordinate \
      INPUT=$input_line \
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
