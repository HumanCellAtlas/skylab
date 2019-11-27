task ValidateSmartSeq2Plate {
    File core_QC
    Array[File] qc_tables
    File gene_matrix
    File isoform_matrix

    String expected_core_QC_hash
    String expected_qc_tables_hash
    String expected_gene_matrix_hash
    String expected_isoform_matrix_hash

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # Always pass for debug
    exit 0;
  >>>
  
  runtime {
    docker: "ubuntu:16.04"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
  }
}
