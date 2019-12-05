task ValidateSmartSeq2Plate {
    Array[File] bam_files
    Array[File] bam_index_files
    Array[File] zarrout

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # Always pass -- just a plumbing test
    exit 0;
  >>>
  
  runtime {
    docker: "ubuntu:16.04"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
  }
}
