task SmartSeq2ZarrConversion {
  # converts some of the outputs of Smart Seq 2 pipeline into a zarr file

  String zarr_output    # output zarr file
  String smartseq_output_folder   # location of the smart seq 2 output
  String sample_id   # id of the sample

  command {
    python3 create_zarr_ss2.py \
      --analysis_output_path ${smartseq_output_folder} \
      --output_path_for_zarr  ${zarr_output} \
      --sample_id ${sample_id} \
      --format DirectoryStore
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.6"
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    File zarr_output = zarr_output
  }
}