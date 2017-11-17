
task SplitBamByCellBarcode {
  File bam_input
  Float size_in_mb = 1024.0

  command {
    SplitBam \
      --bamfile ${bam_input} \
      --output-prefix subfile \
      --subfile-size ${size_in_mb} \
      --tag CB
  }
  
  runtime {
    docker: "humancellatlas/python3-scientific:0.1.3"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 100 HDD"
  }
  
  output {
    Array[File] bam_output_array = glob("subfile_*")
  }
}
