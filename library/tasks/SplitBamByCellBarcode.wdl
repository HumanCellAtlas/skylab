
task SplitBamByCellBarcode {
  File bam_input
  Float size_in_mb = 1024.0
  Int estimated_required_disk = ceil(size(bam_input, "G") * 2.2)

  command {
    SplitBam \
      --bamfile ${bam_input} \
      --output-prefix subfile \
      --subfile-size ${size_in_mb} \
      --tag CB
  }
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.5"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    Array[File] bam_output_array = glob("subfile_*")
  }
}
