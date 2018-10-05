task SmartSeq2ZarrConversion {

  #runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.6_zarr_test"
  File rsem_gene_results # the gene count file "<sample_id>_rsem.genes.results" in the task results folder call-RSEMExpression
  Array[File] smartseq_qc_files # file named "<sample_id>_QCs.csv" in the folder  "call-GroupQCOutputs/glob-*" of the the SS2  output
  String sample_name   # name of the sample

  meta {
    description: "This  task will converts some of the outputs of Smart Seq 2 pipeline into a zarr file"
  }

  parameter_meta {
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    max_retries: "(optional) retry this number of times if task fails -- use with caution, see skylab README for details"
  }

  command {
    set -e

    cp /tools/create_zarr_ss2.py ./

    python3 create_zarr_ss2.py \
       --qc_analysis_output_files_string ${sep=',' smartseq_qc_files} \
       --rsem_genes_results  ${rsem_gene_results} \
       --output_path_for_zarr  "${sample_name}.zarr" \
       --sample_id ${sample_name} \
       --format DirectoryStore


    #get all the files in the zarr folder in  a list
    a=`find "${sample_name}.zarr"  -type f`
    for f in $a; do
       #replace all / to ! as a work around for now.
       b=`echo $f | tr "/" "\!"`
       mv $f $b
    done
  }

  runtime {
    docker: docker
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    Array[File] zarr_output_files = glob("*zarr*")
  }
}
