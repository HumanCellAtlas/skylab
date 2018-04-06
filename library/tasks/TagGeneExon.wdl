task TagGeneExon {
  # this task requires a gtf file with no internal comments
  # each record must have a gene_name and transcript_name in addition to a gene_id and transcript_id
  # the annotation file must be in gtf format (terminal semicolon) not gff format.
  # the gtf must not have any terminal white space at the end of any line
  File annotations_gtf
  File bam_input
  Int estimated_required_disk = ceil((size(bam_input, "G") + size(annotations_gtf, "G")) * 2)

  command {
    set -e

    java -Xms5000m -jar /usr/dropseq-tools-hca/jar/dropseq.jar TagMultiMappedReadWithGeneExon \
      INPUT=${bam_input} \
      OUTPUT=bam_with_gene_exon.bam \
      SUMMARY=gene_exon_tag_summary.log \
      TAG=GE \
      UMI_TAG=UR \
      ALLOW_MULTI_GENE_READS=true \
      FUNCTION_TAG=XF \
      DELETE_SECONDARY_ALIGNMENTS=true \
      ADD_GENE_TAG_TO_INTRONIC_READS=true \
      ANNOTATIONS_FILE=${annotations_gtf}
  }

  # Larger genomes (mouse-human) require a 7.5gb instance; single-organism genomes work with 3.75gb
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/dropseq-tools-hca:1.3.1"
    cpu: 1
    memory: "7.5 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    File bam_output = "bam_with_gene_exon.bam"
    File log = "gene_exon_tag_summary.log"
  }
}
