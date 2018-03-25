task TagGeneExon {
  File annotations_gtf
  File bam_input

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
  # multiply input size by 2
  Int disk = select_first([opt_disk, ceil((size(bam_input, "G") + size(annotations_gtf, "G")) * 2)])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    annotations_gtf: "each record must have a gene_name and transcript_name in addition to a gene_id and transcript_id.  Must be in gtf format"
    bam_input: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
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
