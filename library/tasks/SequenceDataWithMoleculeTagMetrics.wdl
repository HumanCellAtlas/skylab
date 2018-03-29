task CalculateGeneMetrics {
  File bam_input
  Int estimated_required_disk = ceil(size(bam_input, "G") * 2)

  command {
    set -e

    CalculateGeneMetrics -i "${bam_input}" -o gene-metrics.csv.gz
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-sctools:0.1.9"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }

  output {
    File gene_metrics = "gene-metrics.csv.gz"
  }
}


task CalculateCellMetrics {
  File bam_input
  Int estimated_required_disk = ceil(size(bam_input, "G") * 2)

  command {
    set -e

    CalculateCellMetrics -i "${bam_input}" -o cell-metrics.csv.gz
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-sctools:0.1.9"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }

  output {
    File cell_metrics = "cell-metrics.csv.gz"
  }
}


task MergeGeneMetrics {
  Array[File] metric_files

  command {
    set -e

    MergeGeneMetrics -o merged-gene-metrics.csv.gz ${sep=' ' metric_files}
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-sctools:0.1.9"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 20 HDD"
  }

  output {
    File gene_metrics = "merged-gene-metrics.csv.gz"
  }
}


task MergeCellMetrics {
  Array[File] metric_files

  command {
    set -e

    MergeCellMetrics -o merged-cell-metrics.csv.gz ${sep=' ' metric_files}
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-sctools:0.1.9"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 20 HDD"
  }

  output {
    File cell_metrics = "merged-cell-metrics.csv.gz"
  }
}
