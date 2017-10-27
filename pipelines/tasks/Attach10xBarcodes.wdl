
task Attach10xBarcodes {
  File r1  # forward fastq file; contains umi, cell barcode
  File i7  # index fastq file; contains sample barcode
  File u2  # reverse bam file; contains alignable genomic information

  # estimate disk requirements
  Float bam_size = size(u2, "G")
  Float fastq_size = size(r1, "G") + size(i7, "G")
  Int estimated_required_disk = ceil(fastq_size + bam_size * 2.2)

  command {
    Attach10xBarcodes \
      --r1 "${r1}" \
      --i7 "${i7}" \
      --u2 "${u2}" \
      --output-bamfile barcoded.bam
  }
  
  runtime {
    docker: "humancellatlas/python3-scientific:0.1.0"
    cpu: 1
    memory: "3.5 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    File barcoded_bam = "barcoded.bam"
  }
}
