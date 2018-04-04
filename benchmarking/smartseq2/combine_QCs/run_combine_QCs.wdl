task CombineQCs{
  File datafile
  File metfile
  File gtf_file
  Int nthreshold
  String output_name
  String src_dir
  command {
    set -e 
    Rscript /usr/local/scripts/combine_QC_quant_metrics.R  --datafile ${datafile} --metrics ${metfile} --nthreshold ${nthreshold} --gtf ${gtf_file} --out ${output_name}
  }
  output {
    File combined_results = "${output_name}_metrics_combined.csv"
  }
  runtime {
    docker:"gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory:"3.75 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  }
}

workflow RunCombineQC{
  String scripts_dir
  String datafile
  String metfile
  String output_name
  String gtf_file
  Int nthreshold
  call CombineQCs {
    input:
      datafile = datafile,
      metfile = metfile,
      output_name = output_name,
      src_dir = scripts_dir,
      nthreshold = nthreshold,
      gtf_file = gtf_file
  }
  output {
    File CombinedQC = CombineQCs.combined_results
  }
}
