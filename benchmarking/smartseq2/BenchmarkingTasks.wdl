task RunComparativeAnalysis {
  File base_datafile
  File updated_datafile
  File metadata_file
  String metadata_keys
  String output_name

  command {
    set -e
    cp  /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("Compare_data_matrix.R", output_file="${output_name}_Compare_data_matrix.html")' --args \
      "--matrix1=${base_datafile} --matrix2=${updated_datafile} --metaKeys=${metadata_keys} --metadata_file=${metadata_file} --output_prefix=${output_name}"
  }
  output {
    File html = "${output_name}_Compare_data_matrix.html"
    File membership = "${output_name}_tsne_SNN_cluster_data_matrix_comparison.csv"
  }
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  }
}
task RunReproducibilityAnalysis{
  File base_datafile
  File updated_datafile
  File metadata_file
  File gtf_file
  String output_name
  
  command {
    set -e
    cp /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("data_matrix_reproducibility.R", output_file="${output_name}_data_matrix_reproducibility.html")' --args \
      "--matrix1=${base_datafile} --matrix2=${updated_datafile} --metadata_file=${metadata_file} --output_prefix=${output_name} --gtf_file=${gtf_file}"
  }
  output {
    File html = "${output_name}_data_matrix_reproducibility.html"
    File ttest_base = "${output_name}_base_ttest_res.csv"
    File ttest_updated = "${output_name}_updated_ttest_res.csv"
  }
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  } 
}
task RunQCMetricsAnalysis{
  File base_metrics
  File updated_metrics
  String metKeys
  String output_name
  
  command {
    set -e
    cp /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("QC_metrics_analysis.R", output_file="${output_name}_QC_metrics_analysis.html")' --args \
      "--metrics1=${base_metrics} --metrics2=${updated_metrics} --metKeys=${metKeys} --output_prefix=${output_name}"
  }
  output {
    File html = "${output_name}_QC_metrics_analysis.html"
    File test_qc = "${output_name}_tests_stats.csv"
  }
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  } 

}
task RunConfoundingFactorAnalysis{
  File base_metrics
  File updated_metrics
  File metadata_file
  File base_datafile
  File updated_datafile
  String metaKeys
  String output_name
  Int npcs
  
  command {
    set -e
    cp /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("Confounding_factor_analysis.R", output_file="${output_name}_Confounding_factor_analysis.html")' --args \
      "--matrix1=${base_datafile} --matrix2=${updated_datafile} --metadata_file=${metadata_file} --output_prefix=${output_name} --npcs=${npcs} --metaKeys=${metaKeys} --metrics1=${base_metrics} --metrics2=${updated_metrics}"
  } 
  output {
    File html = "${output_name}_Confounding_factor_analysis.html"
    File test_confoundingfactors = "${output_name}_confounding_factors_vars.csv"
  }
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  } 
}

