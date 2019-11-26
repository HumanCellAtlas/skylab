
task AggregateDataMatrix{
  Array[File] filename_array
  String col_name
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-ss2-plate-aggregation:0.0.1"

  command {
    set -e
    python /tools/MergeDataMatrix.py -f ${sep=' ' filename_array}  -t ${col_name} -o ${output_name}
  }

  output{
    File aggregated_result = "${output_name}"
  }

  runtime {
    docker: docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 5
    maxRetries: 1
  }
}

task AggregateQCMetrics{
  Array[File] metric_files
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-ss2-plate-aggregation:0.0.1"
  String run_type
  
  command {
    set -e
    git clone --branch jx-ss2-for-unity https://github.com/HumanCellAtlas/skylab
    python /tools/AggregateMetrics.py -f ${sep=' ' metric_files}   -o ${output_name} -t ${run_type}
  }
  
  output{
    File aggregated_result = output_name+".csv"
  }
  
  runtime {
    docker: docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 5
    maxRetries: 1
  }
}

task AggregateQCMetricsCore{
  Array[File] picard_metric_files
  Array[File] rsem_stats_files
  Array[File] hisat2_stats_files
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-ss2-plate-aggregation:0.0.1"
  String run_type
  
  command {
    set -e
    python /tools/AggregateMetrics.py -f ${sep=' ' picard_metric_files} ${sep=' ' hisat2_stats_files} ${sep=' ' rsem_stats_files}   -o ${output_name} -t ${run_type}
  }
  
  output{
    File aggregated_result = output_name+".csv"
  }
 
  runtime {
    docker: docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 5
    maxRetries: 1
  }
  
}
