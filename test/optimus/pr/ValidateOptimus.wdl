task ValidateOptimus {
      File bam
      File matrix
      File matrix_row_index
      File matrix_col_index
      File cell_metrics
      File gene_metrics

      Int required_disk = ceil((size(bam, "G") + size(matrix, "G")) * 1.1)

      String expected_bam_hash
      String expected_matrix_hash
      String expected_matrix_row_hash
      String expected_matrix_col_hash
      String expected_cell_metric_hash
      String expected_gene_metric_hash

  command <<<

    set -eo pipefail

    ## Generate md5sum for the matrix.npz file, the file is first deflated
    ## to avoid changes in the compression process from affecting output
    ## the process is performed in new directory to avoid other *.npy files
    ## that may be inputs from affecting output
    mkdir matrix_deflate
    cp "${matrix}" matrix_deflate/matrix.npz
    cd matrix_deflate
    unzip matrix.npz
    rm matrix.npz
    matrix_hash=$(find . -name "*.npy" -type f -exec md5sum {} \; | sort -k 2 | md5sum | awk '{print $1}')
    cd ..

    # check matrix row and column indexes files hash
    matrix_row_index_hash=$(zcat "${matrix_row_index}" | md5sum | awk '{print $1}')
    matrix_col_index_hash=$(zcat "${matrix_col_index}" | md5sum | awk '{print $1}')
    
    # check gene and cell metrics
    gene_metric_hash=$(zcat "${gene_metrics}" | md5sum | awk '{print $1}')
    cell_metric_hash=$(zcat "${cell_metrics}" | md5sum | awk '{print $1}')

    # calculate hash as above, but ignore run-specific bam headers
    bam_hash=$(samtools view "${bam}" | md5sum | awk '{print $1}')

    # test each output for equality, echoing any failure states to stdout
    fail=false

    if [ "$matrix_hash" != "${expected_matrix_hash}" ]; then
      >&2 echo "matrix_hash ($matrix_hash) did not match expected hash (${expected_matrix_hash})"
      fail=true
    fi

    if [ "$matrix_row_index_hash" != "${expected_matrix_row_hash)" ]; then
      >&2 echo "matrix_row_index_hash ($matrix_row_index_hash) did not match expected hash (${expected_matrix_row_index_hash})"
      fail=true
    fi

    if [ "$matrix_col_index_hash" != "${expected_matrix_col_index_hash}" ]; then
      >&2 echo "matrix_col_index_hash ($matrix_col_index_hash) did not match expected hash (${expected_matrix_col_index_hash})"
      fail=true
    fi

    if [ "$bam_hash" != "${expected_bam_hash}" ]; then
      >&2 echo "bam_hash ($bam_hash) did not match expected hash (${expected_bam_hash})"
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

    if [ $fail == "true" ]; then exit 1; fi

  >>>
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${required_disk} HDD"
  }
}
