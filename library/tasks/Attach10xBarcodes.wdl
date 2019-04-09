task Attach10xBarcodes {
  File r1_fastq
  File? i1_fastq
  File r2_unmapped_bam
  File whitelist

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-sctools:v0.3.2"
  Int machine_mem_mb = 7500
  Int cpu = 2
  # estimate that bam is approximately the size of all inputs plus 50%
  Int disk = ceil((size(r2_unmapped_bam, "G") + size(r1_fastq, "G") + if (defined(i1_fastq)) then size(i1_fastq, "G") else 0) * 2.5)
  # by default request non preemptible machine to make sure the slow attach barcodes step completes
  Int preemptible = 0

  meta {
    description: "attaches barcodes found in r1 (forward) and i1 (index) fastq files to corresponding reads in the r2 (reverse) bam file"
  }

  parameter_meta {
    r1_fastq: "forward fastq file; contains umi, cell barcode"
    i1_fastq: "optional, index fastq file; contains sample barcode"
    r2_unmapped_bam: "reverse unmapped bam file; contains alignable genomic information"
    whitelist: "10x genomics cell barcode whitelist for 10x V2"
    docker: "(optional) the docker image containing the runtime environment for this task"
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
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
