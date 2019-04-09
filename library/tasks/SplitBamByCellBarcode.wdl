task SplitBamByCellBarcode {
  File bam_input
  Float size_in_mb = 1024.0

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.2"
  Int machine_mem_mb = 3500
  Int cpu = 1
  # estimate that bam is approximately equal in size to the input bam, add 20% buffer
  Int disk = ceil(size(bam_input, "G") * 2.2)
  # by default request non preemptible machine to make sure the slow cell barcode split step completes
  Int preemptible = 0

  meta {
    description: "Splits a bam file into chunks of size_in_mb, guaranteeing that all information for each cell is fully contained in only one of the chunks"
  }

  parameter_meta {
    bam_input: "input bam file"
    size_in_mb: "target size for each output chunk"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
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
