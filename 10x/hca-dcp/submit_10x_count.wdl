task submit {
  String i

  command {
    cat ${i} > out.txt
    echo s >> out.txt
  }
  runtime {
    docker: "ubuntu:17.04"
  }
  output {
    File out = "out.txt"
  }
}

workflow SubmitWrapper10xCount {
  File attach_bcs_and_umis_summary
  File filter_barcodes_summary
  File count_genes_summary
  File extract_reads_summary
  File mark_duplicates_summary
  File raw_gene_bc_matrices_mex
  File raw_gene_bc_matrices_h5
  File filtered_gene_bc_matrices_mex
  File filtered_gene_bc_matrices_h5
  File bam_output
  String i = "asdf"

  call submit {
    input:
      i = i
  }

  output {
    File o = submit.out
  }
}