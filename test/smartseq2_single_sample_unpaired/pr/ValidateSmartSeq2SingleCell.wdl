task ValidateSmartSeq2SingleCell {
      File gene_counts
      String expected_gene_counts_hash

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. We parse the first 7 columns because a bug in RSEM
    # makes later columns non-deterministic.
    gene_counts_hash=$(cut -f 1-7 "${gene_counts}" | md5sum | awk '{print $1}')


    if [ "$counts_hash" != "${expected_gene_counts_hash}" ]; then
      >&2 echo "counts_hash ($counts_hash) did not match expected hash (${expected_counts_hash})"
      fail=true
    fi

    if [ $fail == "true" ]; then exit 1; fi

  >>>
  
  runtime {
    docker: "ubuntu:16.04"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
  }
}
