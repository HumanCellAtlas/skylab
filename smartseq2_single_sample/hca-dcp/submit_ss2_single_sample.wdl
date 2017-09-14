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

workflow SubmitSs2RsemSingleSample {
  File bam_file
  File bam_trans
  File rna_metrics
  File aln_metrics
  File rsem_gene_results
  File rsem_isoform_results
  File rsem_gene_count
  File gene_unique_counts
  File exon_unique_counts
  File transcript_unique_counts
  String i = "asdf"

  call submit {
    input:
      i = i
  }

  output {
    File o = submit.out
  }
}
