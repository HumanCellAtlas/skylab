
task TagGeneExon {
  # this task requires a gtf file with no internal comments
  # each record must have a gene_name and transcript_name in addition to a gene_id and transcript_id
  # the annotation file must be in gtf format (terminal semicolon) not gff format.
  # the gtf must not have any terminal white space at the end of any line
  File annotations_gtf
  File bam_input
  Int estimated_required_disk = ceil((size(bam_input, "G") + size(annotations_gtf, "G")) * 2)

 command {
    TagReadWithGeneExon \
      INPUT=${bam_input} \
      OUTPUT=bam_with_gene_exon.bam \
      SUMMARY=gene_exon_tag_summary.log \
      TAG=GE \
      ANNOTATIONS_FILE=${annotations_gtf}
  }

  # Larger genomes (mouse-human) require a 7.5gb instance; single-organism genomes work with 3.75gb
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-dropseqtools:v0.2.2-1.12"
    cpu: 1
    memory: "7.5 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    File bam_output = "bam_with_gene_exon.bam"
    File log = "gene_exon_tag_summary.log"
  }
} 