version 1.0

task ValidateSmartSeq2Plate {
  input {
    File loom_output
    File truth_loom

    Int disk_size = ceil(size(loom_output,"GiB") + size(truth_loom, "GiB") + 10)
   }

  command <<<

    # catch intermittent failures
    set -eo pipefail

    /tools/loomCompare.py --truth-loom ~{truth_loom} --check-loom ~{loom_output} --delta-cutoff 10

  >>>
  
  runtime {
    docker: "quay.io/humancellatlas/loom-delta-test:0.0.1"
    cpu: 1
    memory: "8 GiB"
    disks: "local-disk 1${disk_size} HDD"
  }
}
