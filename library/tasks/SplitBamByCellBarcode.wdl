task SplitBamByCellBarcode {
  File bam_input
  Float size_in_mb = 1024.0

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.5"
  Int machine_mem_mb = 3000
  Int cpu = 1
  # estimate that bam is approximately equal in size to the input bam, add 20% buffer
  Int disk = ceil(size(bam_input, "G") * 2.2)
  Int preemptible = 0

  meta {
    description: "Splits a bam file into chunks of size_in_mb, guaranteeing that all information for each cell is fully contained in only one of the chunks"
  }

  parameter_meta {
    bam_input: "input bam file"
    size_in_mb: "target size for each output chunk"
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    SplitBam \
      --bamfile ${bam_input} \
      --output-prefix subfile \
      --subfile-size ${size_in_mb} \
      --tags CB CR
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    Array[File] bam_output_array = glob("subfile_*")
  }
}
