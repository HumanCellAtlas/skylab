task SamToFastq {
  File input_bam
  String output_basename

  Int? disk

  command {
    java -Xms2g -jar /root/picard.jar SamToFastq \
    I=${input_bam} \
    F=${output_basename}.fastq
  }
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    # fast q shouldn't be larger than the bam size
    disks: "local-disk " + if defined(disk) then disk + " HDD" else (ceil(size(input_bam, "GB")) * 2) + 10 + " HDD"
    memory: "3 GB"
  }
  output {
    File fastq = "${output_basename}.fastq"
  }
}
