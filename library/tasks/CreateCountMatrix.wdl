
task DropSeqToolsDigitalExpression {
  File bam_input
  File whitelist
  Int estimated_required_disk = ceil((size(bam_input, "G") + size(whitelist, "G")) * 4.0) + 1

  command {
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
    docker: "quay.io/humancellatlas/secondary-analysis-dropseqtools:v0.2.2-1.12"
    cpu: 1
    memory: "7.5 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    File matrix_output = "digital_expression.txt.gz"
    File matrix_summary = "digital_expression_summary.txt"
  }
}
