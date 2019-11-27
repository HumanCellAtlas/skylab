task AggregateDataMatrix {
  Array[File] filename_array
  String col_name
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-ss2-plate-aggregation:0.0.1"

  meta {
    description: "aggregate output count matrix"
  }

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

task AggregateQCMetrics {
  Array[File] metric_files
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-ss2-plate-aggregation:0.0.1"
  String run_type

  meta {
    description: "aggregate count data"
  }
  
  command {
    set -e
    python /tools/AggregateMetrics.py -f ${sep=' ' metric_files}  -o ${output_name} -t ${run_type}
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

task AggregateSmartSeq2Zarr {
    Array[Array[File]?] zarr_input
    String output_file_name
    String docker = "quay.io/humancellatlas/secondary-analysis-ss2-plate-aggregation:0.0.1"
    Int disk = 100

    meta {
      description: "aggregate the zarr output"
    }

    command {
      echo "Hello world" > dummy_output.txt
    }

    output {
        File dummy_output = "dummy_output.txt"
    }

    runtime {
      docker: docker
      memory: "2 GiB"
      disks: "local-disk ${disk} HDD"
      cpu: 1
      preemptible: 3
      maxRetries: 1
    }
}
