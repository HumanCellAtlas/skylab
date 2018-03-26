task CalculateGeneMetrics {
  File bam_input

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:0.1.8"
  Int machine_mem_mb = 3500
  Int cpu = 1
  Int disk = ceil(size(bam_input, "G") * 2)
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

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:0.1.8"
  Int machine_mem_mb = 3500
  Int cpu = 1
  Int disk = ceil(size(bam_input, "G") * 2)
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

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:0.1.8"
  Int machine_mem_mb = 3500
  Int cpu = 1
  Int disk = 20
  Int preemptible = 0

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    metric_files: ""
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
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

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:0.1.8"
  Int machine_mem_mb = 3500
  Int cpu = 1
  Int disk = 20
  Int preemptible = 0

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    metric_files: ""
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
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
