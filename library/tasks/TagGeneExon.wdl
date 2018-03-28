task TagGeneExon {
  File annotations_gtf
  File bam_input

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-dropseqtools:v0.2.2-1.12"
  Int machine_mem_mb = 7500
  Int cpu = 1
  Int disk = ceil((size(bam_input, "G") + size(annotations_gtf, "G")) * 2)
  Int preemptible = 0

  meta {
    description: "Tags any read in bam_input that overlaps an intron or exon interval with the gene that those interals correspond to."
  }

  parameter_meta {
    annotations_gtf: "GTF annotation file for the species that bam input is derived from. Each record must have a gene_name and transcript_name in addition to a gene_id and transcript_id, no white space at the end of any record and must be in gtf format."
    bam_input: "Aligned bam file."
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
  }

 command {
    set -e

    TagReadWithGeneExon \
      INPUT=${bam_input} \
      OUTPUT=bam_with_gene_exon.bam \
      SUMMARY=gene_exon_tag_summary.log \
      TAG=GE \
      ANNOTATIONS_FILE=${annotations_gtf}
  }

  # Larger genomes (mouse-human) require a 7.5gb instance; single-organism genomes work with 3.75gb
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "bam_with_gene_exon.bam"
    File log = "gene_exon_tag_summary.log"
  }
} 
