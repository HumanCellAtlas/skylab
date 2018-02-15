task mergeTPM{
  String target_dir
  String run_name
  String output_name
  String uuid
  String value_name
  command{
    gsutil cp -r ${target_dir} ./   
    python pipeline_testing_scripts/merge_table_quants.py  ${uuid} ${run_name} ${value_name} ${output_name}
  }

  output{
    File merged_tpm_csv = "${output_name}"
  }
  runtime{
    docker:"gcr.io/broad-dsde-mint-dev/analysis-tools:0.0.3"
    memory:"7.5 GB"
    disks: "local-disk 100 HDD"
  }
}

workflow run_merge{
  String scripts_dir
  String uuid
  String value_name
  Array[String] run_names
  Array[String] output_names
  ##Array[String] run_names=['run_pipelines.RunHisat2RsemPipeline.rsem_gene_results','run_pipelines.RunStarPipeline.rsem_gene_results']
  ## Array[String] merged_outputs = ['HISAT2RSEM_gene_TPM.csv','STAR2RSEM_gene_TPM.csv']
  scatter(idx in range(length(run_names))){
    call mergeTPM{
      input:
        target_dir = scripts_dir,
        uuid = uuid,
        run_name = run_names[idx],
        output_name = output_names[idx],
        value_name = value_name
  }
}
}
