
task TagGeneExon {
  File bam_input
  File annotations_gtf
  
  command {
    TagReadWithGeneExon \
      INPUT=${bam_input} \
      OUTPUT=bam_with_gene_exon.bam \
      SUMMARY=gene_exon_tag_summary.log \
      TAG=GE \
      ANNOTATIONS_FILE=${annotations_gtf}
  }
  
  runtime {
    docker: "humancellatlas/dropseqtools:1.12"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 100 HDD"
  }
  
  output {
    File bam_output = "bam_with_gene_exon.bam"
    File log = "gene_exon_tag_summary.log"
  }
} 