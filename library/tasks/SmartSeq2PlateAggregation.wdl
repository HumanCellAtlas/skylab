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
      for f in $( cat input_files.txt )
      do
        mv $f packed_zarr/$( basename $f )
      done

      # Unpack the zarr files by moving them to the appropriate location
      mkdir unpacked_zarr
      /tools/unpackZARR.sh -m -i packed_zarr -o unpacked_zarr

      # Merge the zarr files
      /tools/ss2_zarr_merge.py --input-zarr-directory unpacked_zarr --output-zarr-file output_zarr --plate-sample-id plateid --check-all-headers

      # Flatten the zarr
      mkdir zarrout
      exportfiles=`find output_zarr -type f`
      for f in $exportfiles; do
        newfilename=`echo $f | tr "/" "!"`
        mv $f zarrout/$newfilename
      done

    }

    output {
        Array[File] zarr_output_files = glob("zarrout/*zarr*")
    }

    runtime {
      docker: docker
      memory: "3.5 GiB"
      disks: "local-disk ${disk} HDD"
      cpu: 1
      preemptible: 3
      maxRetries: 1
    }
}
