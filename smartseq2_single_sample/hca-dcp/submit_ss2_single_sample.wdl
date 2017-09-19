task submit {
  String submit_url
  File analysis_json

  command <<<
    submit -submit_url ${submit_url} -analysis_json_path ${analysis_json}
  >>>
  runtime {
    docker: "humancellatlas/secondary-analysis-python:test_submit"
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
  String submit_url
  File analysis_json

  call submit {
    input:
      submit_url = submit_url,
      analysis_json = analysis_json
  }

  output {
    File o = submit.out
  }
}
