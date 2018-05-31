task RunComparativeAnalysis {
  
  meta {
    description: "Run SmartSeq2 Comparative test. In the test, we compare SNN-Cliq clustering results of two data matrix by using Rand adjust index."
  }
  File base_datafile
  File updated_datafile
  File metadata_file
  String metadata_keys
  String output_name

  parameter_meta {
    base_datafile: "data matrix,count  or TPM matrix, of one pipeline"
    updated_datafile: "data matrix, count or TPM matrix, of second pipeline"
    metadata_file: "meta file"
    metadata_keys: "keys in metadata to be used to used as biological labels, such as cell type, cell lineage"
    output_name: "output's prefix"
  }
  
  command {
    set -e
    cp  /usr/local/scripts/*.R ./
    Rscript -e 'rmarkdown::render("Compare_data_matrix.R", output_file="${output_name}_Compare_data_matrix.html")' --args \
      "--matrix1=${base_datafile} --matrix2=${updated_datafile} --metaKeys=${metadata_keys} --metadata_file=${metadata_file} --output_prefix=${output_name}"
  }

  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "15 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  }

  output {
    File html = "${output_name}_Compare_data_matrix.html"
    File membership = "${output_name}_tsne_SNN_cluster_data_matrix_comparison.csv"
  }
}

task RunGeneQuantificationAnalysis{
 
  meta {
    description: "Run SmartSeq2 Gene quantification summary. Pairwise comparison between two data matrix, report the discrepency if gene is low expressed(<low_cut) but high expressed(>high_cut) in anotehr"
  }
 
  File base_datafile
  File updated_datafile
  File gtf_file
  String output_name
  String low_cut
  String high_cut
  
  parameter_meta {
    base_datafile: "data matrix,count  or TPM matrix, of one pipeline"
    updated_datafile: "data matrix, count or TPM matrix, of second pipeline"
    output_name: "output's prefix"
    gtf_file: "gene annotaiton file, in gtf format"
    low_cut: "low threshold cutoff to be used to filter out cells"
    high_cut: "high threshold cutoff to be used to filter out cells"
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
  
  meta {
    description: "Run SmartSeq2 Reproducibility test. In this test, reproducibility between single cell and bulk sample is carried out and the differential genes between single cell and bulk samples will be comparied between to data matrix."
  }

  File base_datafile
  File updated_datafile
  File metadata_file
  File gtf_file
  String output_name
  String groups
  parameter_meta {
    base_datafile: "data matrix,count  or TPM matrix, of one pipeline"
    updated_datafile: "data matrix, count or TPM matrix, of second pipeline"
    output_name: "output's prefix"
    metadata_file: "meta file"
    groups: "labels to be used to identify conditions, such single cell vs bulk, donor1 vs donor2, control vs dieased"
    gtf_file: "gene annotaiton file, in gtf format"
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
 
  meta {
    description: "Run SmartSeq2 QC metrics Comparison analysis. 3 statistical tests are carried out to compare QC metrics: a Wilcoxon Pairwise test, a linear Regression model test and a  Kolmogorov–Smirnov test "
  }
 
  File base_metrics
  File updated_metrics
  String met_keys
  String output_name
  
  parameter_meta {
    base_metrics: "QC metrics, each row represents a metric and each column represents a cell"
    updated_metrics: "QC metrics, each row represents a metric and each column represents a cell"
    output_name: "output's prefix"
    met_keys: "names of QC metrics to extract and compare"
  }
  
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
  
  meta {
    description: This analysis is to evaluate and compare confounding factors impacts of two pipelines.
  }

  File base_metrics
  File updated_metrics
  File metadata_file
  File base_datafile
  File updated_datafile
  String meta_keys
  String output_name
  Int npcs
  
  parameter_meta {
    base_metrics: "QC metrics, each row represents a metric and each column represents a cell"
    updated_metrics: "QC metrics, each row represents a metric and each column represents a cell"
    base_datafile: "data matrix,count  or TPM matrix, of one pipeline"
    updated_datafile: "data matrix, count or TPM matrix, of second pipeline"
    metadata_file: "meta file"
    output_name: "output's prefix"
    meta_keys: "names of meta values to extract and use as biological confounding factors"
    npcs: "number of PCs to extract"
  } 

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
