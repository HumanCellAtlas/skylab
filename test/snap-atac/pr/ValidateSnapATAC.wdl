version 1.0

task ValidateSnapATAC {
  input {
      File snap
      File snapqc
      File bam

      Int required_disk = ceil((size(bam, "G") + size(snap, "G")) * 1.2)

      String expected_snap_hash
      String expected_snapqc_hash
      String expected_bam_hash
   }

  command <<<

    exit 0;

  >>>
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${required_disk} HDD"
  }
}
