task DataMatrixReproducibility {
  File base_datafile
  File updated_datafile
  File gtf_file
  File metadata_file
  String output_name
  String src_dir
  
  command {
    set -e
    mkdir ${output_name}
    Rscript /usr/local/scripts/data_matrix_reproducibility.R  \
      --matrix1 ${base_datafile} \
      --matrix2 ${updated_datafile} \
      --gtf_file ${gtf_file} \
      --metadata_file ${metadata_file} \
      --output_prefix ${output_name}/${output_name} 
    tar -zcvf ${output_name}_matrix.tar.gz ${output_name}
  }
  output {
    File results = "${output_name}_matrix.tar.gz"
  }
  runtime {
    docker:"gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory:"7.5 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  }
}

workflow RunDataMatrixReproducibility {
  String scripts_dir
  String base_datafile
  String updated_datafile
  String output_name
  String gtf_file
  String metadata_file
  call DataMatrixReproducibility {
    input:
      base_datafile = base_datafile,
      updated_datafile = updated_datafile,
      output_name = output_name,
      gtf_file = gtf_file,
      src_dir = scripts_dir,
      metadata_file = metadata_file
  }
  output {
    File AnalysisResults = DataMatrixReproducibility.results
  }
}
