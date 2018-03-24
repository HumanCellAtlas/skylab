task CellSortBam {
  File bam_input
  Int estimated_required_disk = ceil(size(bam_input, "G") * 8)
  
  command {
    set -e

    samtools sort -t GE -t UB -t CB -o cell-sorted.bam "${bam_input}"
  }
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 2
    memory: "7.5 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    File bam_output = "cell-sorted.bam"
  }
}


task GeneSortBam {
  File bam_input
  Int estimated_required_disk = ceil(size(bam_input, "G") * 4)

  command {
    set -e

    samtools sort -t UB -t CB -t GE -o gene-sorted.bam "${bam_input}"
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 2
    memory: "7.5 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }

  output {
    File bam_output = "gene-sorted.bam"
  }
}
