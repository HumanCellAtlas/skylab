task ValidateSmartSeq2SingleCell {
      File counts
      String expected_counts_hash

      File? target_metrics
      String expected_metrics_hash

      Array[File] fastqc_htmls
      Array[String] expected_fastqc_html_strings

      Int expected_n_fastqc_zips
      Int n_fastqc_zips

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. We parse the first 7 columns because a bug in RSEM
    # makes later columns non-deterministic.
    counts_hash=$(cut -f 1-7 "${counts}" | md5sum | awk '{print $1}')

    # this parses the picard metrics file with awk to remove all the run-specific comment lines (#)
    target_metrics_hash=$(cat "${target_metrics}" | awk 'NF && $1!~/^#/' | md5sum | awk '{print $1}')

    if [ "$counts_hash" != "${expected_counts_hash}" ]; then
      >&2 echo "counts_hash ($counts_hash) did not match expected hash (${expected_counts_hash})"
      fail=true
    fi

    if [ "$target_metrics_hash" != "${expected_metrics_hash}" ]; then
      >&2 echo "target_metrics_hash ($target_metrics_hash) did not match expected hash (${expected_metrics_hash})"
      fail=true
    fi

    #search for expected strings in fastqc html
    for string in "${sep='" "' expected_fastqc_html_strings}"; do
        fail_html=true
        for htmlfile in ${sep=' ' fastqc_htmls}; do
          if grep "$string" $htmlfile; then
            fail_html=false
          fi
        done
        if [ $fail_html == "true" ]; then
          >&2 echo "expected string ($string) not found in fastqc html files"
          fail=true
        fi
    done

    if [ ${expected_n_fastqc_zips} != ${n_fastqc_zips} ]; then
      >&2 echo "number of fastqc zip (${n_fastqc_zips}) did not match expected number (${expected_n_fastqc_zips})"
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
