task Attach10xBarcodes {
  File r1_fastq
  File? i1_fastq
  File r2_unmapped_bam
  File whitelist

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-sctools:0.1.6"])
  Int machine_mem_mb = select_first([opt_memory_gb, 7]) * 1000
  Int cpu = select_first([opt_cpu, 2])
  # estimate that bam is approximately the size of all inputs plus 50%
  Int disk = select_first([opt_disk, ceil((size(r2_unmapped_bam, "G") + size(r1_fastq, "G") + if (defined(i1_fastq)) then size(i1_fastq, "G") else 0) * 2.5)])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "attaches barcodes found in r1 and i1"
  }

  parameter_meta {
    r1_fastq: "forward fastq file; contains umi, cell barcode"
    i1_fastq: "optional, index fastq file; contains sample barcode"
    r2_unmapped_bam: "reverse bam file; contains alignable genomic information"
    whitelist: "10x genomics cell barcode whitelist for 10x V2"
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    Attach10xBarcodes \
      --r1 "${r1_fastq}" \
      ${"--i1 " + i1_fastq} \
      --u2 "${r2_unmapped_bam}" \
      --output-bamfile barcoded.bam \
      --whitelist "${whitelist}"
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "barcoded.bam"
  }
}
