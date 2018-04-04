task AnalysisQCMetrics{
  File base_qc_metrics
  File updated_qc_metrics
  String output_name
  String metrics_keys
  String src_dir 
  command {
    set -e 
    mkdir ${output_name}
    Rscript /usr/local/scripts/QC_metrics_analysis.R --bmetrics ${base_qc_metrics} --umetrics ${updated_qc_metrics} --out ${output_name}/${output_name} --metKeys ${metrics_keys}
    tar -zcvf ${output_name}_picard.tar.gz ${output_name}
  }
  output {
    File results = "${output_name}_picard.tar.gz"
  }
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "3.75 GB"
    disks:  "local-disk 50 HDD"
    preemptible: 5
  }
}

workflow RunQCAnalysis {
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
  output {
    File AnalysisResults = AnalysisQCMetrics.results
  }
}
