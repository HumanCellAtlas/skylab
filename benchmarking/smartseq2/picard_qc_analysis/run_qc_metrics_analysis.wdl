task analysis_QC_metrics{
  File qc_metrics_file1
  File qc_metrics_file2
  String output_name
  String src_dir 
  command {
    set -e 
    mkdir /tmp
    cd /tmp
    gsutil cp -r ${src_dir} ./
    Rscript pipeline_testing_scripts/QC_metrics_analysis.R ${qc_metrics_file1} ${qc_metrics_file2} ${output_name}
    
  }
}

