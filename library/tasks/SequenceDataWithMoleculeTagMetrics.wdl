task CalculateGeneMetrics {
  File bam_input

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.3"
  Int machine_mem_mb = 20000
  Int cpu = 1
  Int disk = ceil(size(bam_input, "G") * 4)
  Int preemptible = 3

  meta {
    description: "Calculate gene metrics from the reads in bam_input."
  }

  parameter_meta {
    bam_input: "An aligned bam file augmented with CB, UB, GE, CY, UY, and XF tags."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
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
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.3"
  Int machine_mem_mb = 3500
  Int cpu = 1
  Int disk = ceil(size(bam_input, "G") * 2)
  Int preemptible = 3

  meta {
    description: "Calculate cell metrics from the reads in bam_input."
  }

  parameter_meta {
    bam_input: "An aligned bam file augmented with CB, UB, GE, CY, UY, and XF tags."
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
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
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.3"
  Int machine_mem_mb = 3500
  Int cpu = 1
  Int disk = 20
  Int preemptible = 3

  meta {
    description: "Merge an array of gene metric files with the same metric categories and potentially overlapping sets of gene features"
  }

  parameter_meta {
    metric_files: "A set of metrics files, each measuring a potentially overlapping set of genes in disjoint sets of cells"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
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
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.3"
  Int machine_mem_mb = 3500
  Int cpu = 1
  Int disk = 20
  Int preemptible = 3

  meta {
    description: "Concatenate multiple cell metrics files into a single matrix"
  }

  parameter_meta {
    metric_files: "An array of cell metrics files that contain the same metric types, but different sets of cells"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
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
