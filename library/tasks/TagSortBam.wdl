task CellSortBam {
  File bam_input

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
  Int machine_mem_mb = 7500
  Int cpu = 2
  Int disk = ceil(size(bam_input, "G") * 8)
  Int preemptible = 0

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_input: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }
  
  command {
    set -e

    samtools sort -t GE -t UB -t CB -o cell-sorted.bam "${bam_input}"
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "cell-sorted.bam"
  }
}

task GeneSortBam {
  File bam_input

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
  Int machine_mem_mb = 7500
  Int cpu = 2
  Int disk = ceil(size(bam_input, "G") * 4)
  Int preemptible = 0

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_input: ""
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    samtools sort -t UB -t CB -t GE -o gene-sorted.bam "${bam_input}"
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "gene-sorted.bam"
  }
}
