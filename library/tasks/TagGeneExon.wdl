
task TagGeneExon {
  File bam_input
  File annotations_gtf
  Int estimated_required_disk = ceil((size(bam_input, "G") + size(annotations_gtf, "G")) * 2)
  
  command {
    TagReadWithGeneExon \
      INPUT=${bam_input} \
      OUTPUT=bam_with_gene_exon.bam \
      SUMMARY=gene_exon_tag_summary.log \
      TAG=GE \
      ANNOTATIONS_FILE=${annotations_gtf}
  }
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-dropseqtools:v0.2.2-1.12"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    File bam_output = "bam_with_gene_exon.bam"
    File log = "gene_exon_tag_summary.log"
  }
} 