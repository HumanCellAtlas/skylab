task DropSeqToolsDigitalExpression {
  File bam_input
  File whitelist
  Int estimated_required_disk = ceil((size(bam_input, "G") + size(whitelist, "G")) * 4.0) + 1

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-dropseqtools:v0.2.2-1.12"
  Int machine_mem_mb = 7500
  Int cpu = 1
  Int disk = ceil((size(bam_input, "G") + size(whitelist, "G")) * 4.0) + 1
  Int preemptible = 0

  meta {
    description: "Constructs a tab-delimited gzipped count matrix from a bam file containing reads marked with cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
  }

  parameter_meta {
    bam_input: "input bam file marked with cell barcodes (CB), molecule barcodes (UB) and gene ids (GE)"
    whitelist: "10x genomics cell barcode whitelist for 10x V2. Only CB found in this list are included in the count matrix"
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    DigitalExpression \
      SUMMARY=digital_expression_summary.txt \
      OUTPUT=digital_expression.txt.gz \
      CELL_BARCODE_TAG=CB \
      MOLECULAR_BARCODE_TAG=UB \
      GENE_EXON_TAG=GE \
      CELL_BC_FILE=${whitelist} \
      USE_STRAND_INFO=false \
      INPUT=${bam_input}
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File matrix_output = "digital_expression.txt.gz"
    File matrix_summary = "digital_expression_summary.txt"
  }
}
