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
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_input: ""
    whitelist: ""
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
