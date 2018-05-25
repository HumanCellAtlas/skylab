task RunComparativeAnalysis {
  File base_datafile
  File updated_datafile
  File metadata_file
  String metadata_keys
  String output_name

  meta {
    description: "Run SmartSeq2 Comparative test. In the test, two data matrix will be tested by their clustering results."
  } 

  command {
    set -e
    cp  /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("Compare_data_matrix.R", output_file="${output_name}_Compare_data_matrix.html")' --args \
      "--matrix1=${base_datafile} --matrix2=${updated_datafile} --metaKeys=${metadata_keys} --metadata_file=${metadata_file} --output_prefix=${output_name}"
  }

  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  }

  output {
    File html = "${output_name}_Compare_data_matrix.html"
    File membership = "${output_name}_tsne_SNN_cluster_data_matrix_comparison.csv"
  }
}

task RunGeneQuantificationAnalysis{
  File base_datafile
  File updated_datafile
  File gtf_file
  String output_name
  String low_cut
  String high_cut
  
  meta {
    description: "Run SmartSeq2 Gene quantification summary. Pairwise comparison between two data matrix, report the discrepency if gene is low expressed(<low_cut) but high expressed(>high_cut) in anotehr"
  } 
  
  command {
    set -e
    cp /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("Compare_Quantification.R",output_file="${output_name}_Compare_Quantification.html")' --args \
      "--matrix1=${base_datafile} --matrix2=${updated_datafile} --output_prefix=${output_name} --gtf_file=${gtf_file} --low_cut=${low_cut} --high_cut=${high_cut}" 
  }   
  
  output {
    File html = "${output_name}_Compare_Quantification.html"
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
  String groups

  meta {
    description: "Run SmartSeq2 Reproducibility test. In this test, reproducibility between single cell and bulk sample is carried out and the differential genes between single cell and bulk samples will be comparied between to data matrix."
  } 

  command {
    set -e
    cp /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("data_matrix_reproducibility.R", output_file="${output_name}_data_matrix_reproducibility.html")' --args \
      "--matrix1=${base_datafile} --matrix2=${updated_datafile} --metadata_file=${metadata_file} --output_prefix=${output_name} --gtf_file=${gtf_file} --groups=${groups}"
  }

  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  } 

  output {
    File html = "${output_name}_data_matrix_reproducibility.html"
    File ttest_base = "${output_name}_base_ttest_res.csv"
    File ttest_updated = "${output_name}_updated_ttest_res.csv"
  }
}

task RunQCMetricsAnalysis{
  File base_metrics
  File updated_metrics
  String met_keys
  String output_name
  
  command {
    set -e
    cp /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("QC_metrics_analysis.R", output_file="${output_name}_QC_metrics_analysis.html")' --args \
      "--metrics1=${base_metrics} --metrics2=${updated_metrics} --metKeys=${met_keys} --output_prefix=${output_name}"
  }

  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  } 

  output {
    File html = "${output_name}_QC_metrics_analysis.html"
    File test_qc = "${output_name}_tests_stats.csv"
  }
}

task RunConfoundingFactorAnalysis{
  File base_metrics
  File updated_metrics
  File metadata_file
  File base_datafile
  File updated_datafile
  String meta_keys
  String output_name
  Int npcs
  
  command {
    set -e
    cp /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("Confounding_factor_analysis.R", output_file="${output_name}_Confounding_factor_analysis.html")' --args \
      "--matrix1=${base_datafile} --matrix2=${updated_datafile} --metadata_file=${metadata_file} --output_prefix=${output_name} --npcs=${npcs} --metaKeys=${meta_keys} --metrics1=${base_metrics} --metrics2=${updated_metrics}"
  } 

  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  } 

  output {
    File html = "${output_name}_Confounding_factor_analysis.html"
    File test_confoundingfactors = "${output_name}_confounding_factors_vars.csv"
  }
}