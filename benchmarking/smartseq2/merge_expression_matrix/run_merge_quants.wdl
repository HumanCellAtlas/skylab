task MergeTPM {
  String target_dir
  String file_name
  String output_name
  String uuid
  String value_name
  command {
    gsutil cp  ${target_dir}/merge_table_quants.py ./   
    python merge_table_quants.py  -u ${uuid} -fn ${file_name} -t ${value_name} -o ${output_name}
  }
  output {
    File merged_tpm_csv = "${output_name}"
  }
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/analysis-tools:0.0.5"
    memory: "7.5 GB"
    disks:  "local-disk 100 HDD"
    preemptible: 5   
  }
}

workflow run_merge {
  String scripts_dir
  String uuid
  String value_name
  Array[String] file_names
  Array[String] output_names
  scatter(idx in range(length(run_names))) {
    call MergeTPM {
      input:
        target_dir = scripts_dir,
        uuid = uuid,
        file_name = file_names[idx],
        output_name = output_names[idx],
        value_name = value_name
  }
}
}
