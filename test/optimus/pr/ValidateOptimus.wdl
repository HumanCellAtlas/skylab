
task ValidateOptimus {
      File bam
      File matrix
      File matrix_summary
      File picard_metrics
      Int required_disk = ceil((size(bam, "G") + size(matrix, "G") + size(matrix_summary, "G") + size(picard_metrics, "G")) * 1.1)

      String expected_bam_hash
      String expected_matrix_hash
      String expected_matrix_summary_hash
      String expected_picard_metrics_hash

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed
    matrix_hash=$(md5sum "${matrix}" | awk '{print $1}')
    matrix_summary_hash=$(md5sum "${matrix_summary}" | awk '{print $1}')

    # calculate hash as above, but ignore run-specific bam headers
    bam_hash=$(samtools view "${bam}" | md5sum | awk '{print $1}')

    # the picard metrics are contained in a .tar.gz file; in addition to the above processing,
    # this unzips the archive, and parses it with awk to remove all the run-specific comment lines (#)
    picard_metrics_hash=$(gunzip -c "${picard_metrics}" | awk 'NF && $1!~/^#/' | md5sum | awk '{print $1}')

    # test each output for equivalence, echoing any failure states to stdout
    fail=false
    if [ "$bam_hash" != "${expected_bam_hash}" ]; then
      >&2 echo "bam_hash ($bam_hash) did not match expected hash (${expected_bam_hash})"
      fail=true
    fi
    
    if [ "$matrix_hash" != "${expected_matrix_hash}" ]; then
      >&2 echo "matrix_hash ($matrix_hash) did not match expected hash (${expected_matrix_hash})"
      fail=true
    fi

    if [ "$matrix_summary_hash" != "${expected_matrix_summary_hash}" ]; then
      >&2 echo "matrix_summary_hash ($matrix_summary_hash) did not match expected hash (${expected_matrix_summary_hash})"
      fail=true
    fi

    if [ "$picard_metrics_hash" != "${expected_picard_metrics_hash}" ]; then
      >&2 echo "picard_metrics_hash ($picard_metrics_hash) did not match expected hash (${expected_picard_metrics_hash})"
      fail=true
    fi
    
    if [ $fail == "true" ]; then exit 1; fi

  >>>
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${required_disk} HDD"
  }
}
