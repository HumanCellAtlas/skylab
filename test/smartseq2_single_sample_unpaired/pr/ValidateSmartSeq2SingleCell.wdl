task ValidateSmartSeq2SingleCell {
      File gene_counts
      String expected_gene_counts_hash

      String dollar='$'

  command <<<

    # catch intermittent failures
    set -xeo pipefail

    # The test data is a cd4+ t cell, so make sure we get a little cd4, cd3d
    # and cd3g
    gene_ids=(ENSG00000160654 ENSG00000167286 ENSG00000010610)
    for gene_id in "${dollar}{gene_ids[@]}"; do
        tpm=$(grep "$gene_id" "${gene_counts}" | cut -f6)
        echo $tpm
        if (( $(echo "$tpm == 0.0" | bc -l) )); then
            >&2 echo "Count not find gene id $gene_id"
            fail=true
        fi
    done
    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. We parse the first 7 columns because a bug in RSEM
    # makes later columns non-deterministic.
    gene_counts_hash=$(cut -f 1-7 "${gene_counts}" | md5sum | awk '{print $1}')

    if [ "$gene_counts_hash" != "${expected_gene_counts_hash}" ]; then
      >&2 echo "counts_hash ($gene_counts_hash) did not match expected hash (${expected_gene_counts_hash})"
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
