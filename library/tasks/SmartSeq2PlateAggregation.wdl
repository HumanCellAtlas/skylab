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

    Array[File] zarrs_flat = flatten(select_all(zarr_input))

    meta {
      description: "aggregate the zarr output"
    }

    command {
      set -e

      # Get all the files we are operating on in a file
      mv ${write_lines(zarrs_flat)} input_files.txt

      # Move all the packed zarrs into a single directory for unpackZarr to work on
      mkdir packed_zarr
      for f in $( cat input_files )
      do
        mv $f packed_zarr/$( basename $f )
      done

      # Unpack the zarr files by moving them to the appropriate location
      mkdir unpacked_zarr
      /tools/unpackZARR.sh -m -i packed_zarr -o unpacked_zarr

      # Check that the directory structure has been unpacked correctly
      tree > dummy_output.txt

      # TODO process the directory structure and generate a single ZARR with all the information
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
