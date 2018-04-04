task MergePicardMetrics {
  String target_dir
  String met_name
  String uuid
  String output_name
  command {
    set -e
    python /usr/local/scripts/merge_picard_mets.py  -u ${uuid} -m ${met_name} -o ${output_name}
  }
  output {
    File merged_metrics = "${output_name}"
  }
  runtime {
    docker:"gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory:"3.75 GB"
    disks: "local-disk 50 HDD"
    preemptible: 5
  }
}
task AppendMetricsFiles { 
  Array[File] met_files
  String output_name
  command {
    set -e
    isheader=0
    for f in ${sep=' ' met_files}; do
      if [[ $isheader == 0 ]]; then
        isheader=1
        cat $f >${output_name}
      else
        tail -n+2 $f >>${output_name}
      fi
    done  
  }
  output {
    File merged_metrics = "${output_name}"
  }
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
    memory: "3.75 GB"
    disks:  "local-disk 10 HDD"
    preemptible: 5
  }
}
workflow RunMergePicardQC {
  String scripts_dir
  String uuid
  Array[String] metrics_names
  String combined_file
  scatter(idx in range(length(metrics_names))) {
    call MergePicardMetrics {
      input:
        target_dir = scripts_dir,
        uuid = uuid,
        met_name = metrics_names[idx],
        output_name = metrics_names[idx],
    }
  }
 call AppendMetricsFiles {
  input:
    met_files = MergePicardMetrics.merged_metrics,
    output_name = combined_file
  }
}
