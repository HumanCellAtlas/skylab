task ValidateOptimus {
      File bam
      File matrix
      File gene_metrics
      File cell_metrics
      Array[File] fastqc_htmls
      Int n_fastqc_zips
      Int required_disk = ceil((size(bam, "G") + size(matrix, "G")) * 1.1)

      String expected_bam_hash
      String expected_matrix_hash
      String expected_gene_metric_hash
      String expected_cell_metric_hash
      Int expected_n_fastqc_zips
      Array[String] expected_fastqc_html_strings

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. gzipped files are unzipped to avoid hashing malleable
    # metadata


    unzip "${matrix}"
    matrix_hash=$(find . -name "*.npy" -type f -exec md5sum {} \; | sort -k 2 | md5sum | awk '{print $1}')
    gene_metric_hash=$(zcat "${gene_metrics}" | md5sum | awk '{print $1}')
    cell_metric_hash=$(zcat "${cell_metrics}" | md5sum | awk '{print $1}')

    # calculate hash as above, but ignore run-specific bam headers
    bam_hash=$(samtools view "${bam}" | md5sum | awk '{print $1}')

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

    if [ "$gene_metric_hash" != "${expected_gene_metric_hash}" ]; then
      >&2 echo "gene_metric_hash ($gene_metric_hash) did not match expected hash (${expected_gene_metric_hash})"
      fail=true
    fi

    if [ "$cell_metric_hash" != "${expected_cell_metric_hash}" ]; then
      >&2 echo "cell_metric_hash ($cell_metric_hash) did not match expected hash (${expected_cell_metric_hash})"
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
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${required_disk} HDD"
  }
}
