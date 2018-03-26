task CalculateGeneMetrics {
  File bam_input

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-sctools:0.1.8"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  Int cpu = select_first([opt_cpu, 1])
  Int disk = select_first([opt_disk, ceil(size(bam_input, "G") * 2)])
  Int preemptible = select_first([opt_preemptible, 0])

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

    CalculateGeneMetrics -i "${bam_input}" -o gene-metrics.csv.gz
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File gene_metrics = "gene-metrics.csv.gz"
  }
}

task CalculateCellMetrics {
  File bam_input

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-sctools:0.1.8"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  Int cpu = select_first([opt_cpu, 1])
  Int disk = select_first([opt_disk, ceil(size(bam_input, "G") * 2)])
  Int preemptible = select_first([opt_preemptible, 0])

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

    CalculateCellMetrics -i "${bam_input}" -o cell-metrics.csv.gz
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File cell_metrics = "cell-metrics.csv.gz"
  }
}


task MergeGeneMetrics {
  Array[File] metric_files

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-sctools:0.1.8"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  Int cpu = select_first([opt_cpu, 1])
  Int disk = select_first([opt_disk, 20])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    metric_files: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    MergeGeneMetrics -o merged-gene-metrics.csv.gz ${sep=' ' metric_files}
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File gene_metrics = "merged-gene-metrics.csv.gz"
  }
}

task MergeCellMetrics {
  Array[File] metric_files

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-sctools:0.1.8"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  Int cpu = select_first([opt_cpu, 1])
  Int disk = select_first([opt_disk, 20])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    metric_files: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    MergeCellMetrics -o merged-cell-metrics.csv.gz ${sep=' ' metric_files}
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File cell_metrics = "merged-cell-metrics.csv.gz"
  }
}
