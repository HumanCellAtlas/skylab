task SmartSeq2ZarrConversion {
  # converts some of the outputs of Smart Seq 2 pipeline into a zarr file

  File rsem_gene_results # the gene count file
  Array[File] smartseq_qc_files   # location of the smart seq 2 output
  String sample_name   # id of the sample

  command {

    set -e

    cp /tools/create_zarr_ss2.py ./

    python3 ./create_zarr_ss2.py  -h

    python3 create_zarr_ss2.py \
       --qc_analysis_output_files_string ${sep=',' smartseq_qc_files} \
       --rsem_genes_results  ${rsem_gene_results} \
       --output_path_for_zarr  "${sample_name}.zip" \
       --sample_id ${sample_name} \
       --format ZipStore
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.6_zarr_test"
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    Array[File]  zarr_output_files = glob("*.zip")
  }
}
