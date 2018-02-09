
task Attach10xBarcodes {
  # attaches barcodes found in r1 and i1

  File r1  # forward fastq file; contains umi, cell barcode
  File u2  # reverse bam file; contains alignable genomic information
  File? i1  # optional, index fastq file; contains sample barcode
  File whitelist  # 10x genomics cell barcode whitelist for 10x V2

  # estimate disk requirements
  Float bam_size = size(u2, "G")
  Float fastq_size = size(r1, "G")  + if (defined(i1)) then size(i1, "G") else 0
  Int estimated_required_disk = ceil((fastq_size + bam_size) * 2.5)

  command {
    Attach10xBarcodes \
      --r1 "${r1}" \
      ${"--i1 " + i1} \
      --u2 "${u2}" \
      --output-bamfile barcoded.bam \
      --whitelist "${whitelist}"
  }
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-sctools:0.1.6"
    cpu: 2
    memory: "7.5 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    File bam_output = "barcoded.bam"
  }
}
