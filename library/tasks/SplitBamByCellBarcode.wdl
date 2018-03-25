task SplitBamByCellBarcode {
  File bam_input
  Float size_in_mb = 1024.0

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.5"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  Int cpu = select_first([opt_cpu, 1])
  # estimate that bam is approximately equal in size to the input bam, add 20% buffer
  Int disk = select_first([opt_disk, ceil(size(bam_input, "G") * 2.2)])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_input: "input bam file"
    size_in_mb: "threshold of how big each bam split should be"
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
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
