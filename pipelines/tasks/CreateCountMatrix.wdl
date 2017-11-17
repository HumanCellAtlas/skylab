# AUTHOR: 
# Created on 2017/11/22

task DropSeqToolsDigitalExpression {
  File bam_input
  File whitelist

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
    docker: "humancellatlas/dropseqtools:1.12"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 100 HDD"
  }
  
  output {
    File matrix_output = "digital_expression.txt.gz"
    File matrix_summary = "digital_expression_summary.txt"
  }
}
