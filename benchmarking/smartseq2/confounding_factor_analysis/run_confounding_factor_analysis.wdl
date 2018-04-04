task AnalysisConfoundingFactors{
  File base_datafile
  File base_metrics
  File updated_datafile
  File updated_metrics
  Int  num_pcs
  String output_name
  String src_dir
  
  command {
    set -e
    mkdir ${output_name}
    Rscript  /usr/local/scripts/Confounding_factor_analysis.R \
      --bdatafile ${base_datafile} \
      --bmetrics ${base_metrics} \
      --udatafile ${updated_datafile} \
      --umetrics ${updated_metrics} \
      --npcs ${num_pcs} \
      --out ${output_name}/${output_name} 
    tar -zcvf ${output_name}_confounding.tar.gz ${output_name}
  }
  
  output {
    File results = "${output_name}_confounding.tar.gz"
  }
  
  runtime {
    docker:"gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory:"7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  }
}

workflow RunConfoundingFactorsAnalysis{
  String scripts_dir
  String base_datafile
  String updated_datafile
  String base_metrics
  String updated_metrics
  String output_name
  Int    npcs

  call AnalysisConfoundingFactors {
    input:
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      base_metrics = base_metrics,
      updated_metrics = updated_metrics,
      output_name = output_name,
      src_dir = scripts_dir,
      num_pcs = npcs,
  }
  output {
    File AnalysisResults = AnalysisConfoundingFactors.results
  }
}
