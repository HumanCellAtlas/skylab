task DropSeqToolsDigitalExpression {
  File bam_input
  File whitelist
  Int estimated_required_disk = ceil((size(bam_input, "G") + size(whitelist, "G")) * 4.0) + 1

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-dropseqtools:v0.2.2-1.12"])
  Int machine_mem_mb = select_first([opt_memory_gb, 7]) * 1000
  Int cpu = select_first([opt_cpu, 1])
  Int disk = select_first([opt_disk, ceil((size(bam_input, "G") + size(whitelist, "G")) * 4.0) + 1])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_input: ""
    whitelist: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
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
