task AnalysisQCMetrics{
  File base_qc_metrics
  File updated_qc_metrics
  String output_name
  String metrics_keys
  String src_dir 
  command {
    set -e 
    mkdir ${output_name}
    gsutil cp ${src_dir}/QC_metrics_analysis.R ./
    Rscript QC_metrics_analysis.R ${base_qc_metrics} ${updated_qc_metrics} ${output_name} ${metrics_keys}
    tar -zcvf ${output_name}.tar.gz ${output_name}
  }
  output {
    File analysis_results = "${output_name}.tar.gz"
  }
  runtime {
    docker:"gcr.io/broad-dsde-mint-dev/analysis-tools:0.0.4"
    memory:"3.75 GB"
    disks: "local-disk 50 HDD"
  }
}
workflow run_analysis {
  String scripts_dir
  String base_qc
  String updated_qc
  String output_dir
  String metrics_keys
  call AnalysisQCMetrics {
    input:
      base_qc_metrics = base_qc,
      updated_qc_metrics = updated_qc,
      output_name = output_dir,
      src_dir = scripts_dir,
      metrics_keys = metrics_keys
  }
}

