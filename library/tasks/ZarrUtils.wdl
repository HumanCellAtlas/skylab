task SmartSeq2ZarrConversion {

  #runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.9"
  # the gene count file "<sample_id>_rsem.genes.results" in the task results folder call-RSEMExpression
  File rsem_gene_results
  # file named "<sample_id>_QCs.csv" in the folder  "call-GroupQCOutputs/glob-*" of the the SS2  output
  Array[File] smartseq_qc_files
  # name of the sample
  String sample_name

  Int preemptible = 3

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
    set -euo pipefail

    python3 /tools/create_zarr_ss2.py \
       --qc_files ${sep=' ' smartseq_qc_files} \
       --rsem_genes_results  ${rsem_gene_results} \
       --output_zarr_path  "${sample_name}.zarr" \
       --sample_id ${sample_name} \
       --format DirectoryStore

    mkdir zarrout
    # get all the files in the zarr folder in  a list
    a=`find "${sample_name}.zarr"  -type f`
    for f in $a; do
       # replace all / to ! as a work around for now.
       b=`echo $f | tr "/" "\!"`
       mv $f zarrout/$b
    done
  }

  runtime {
    docker: docker
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    preemptible: preemptible
  }

  output {
    Array[File] zarr_output_files = glob("zarrout/*zarr*")
  }
}


task OptimusZarrConversion {
  #runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.9"
  # name of the sample
  String sample_id
  # the file "merged-cell-metrics.csv.gz" that contains the cellwise metrics
  File cell_metrics
  # the file "merged-gene-metrics.csv.gz" that contains the  genwise metrics
  File gene_metrics
  # file (.npz)  that contains the count matrix
  File sparse_count_matrix
  # file (.npy) that contains the array of cell barcodes
  File cell_id
  # file (.npy) that contains the array of gene names
  File gene_id
  # emptydrops output metadata
  File empty_drops_result

  Int preemptible = 3

  meta {
    description: "This task will converts some of the outputs of Optimus pipeline into a zarr file"
  }

  parameter_meta {
    machine_mem_mb: "(optional) the amount of memory (MB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    max_retries: "(optional) retry this number of times if task fails -- use with caution, see skylab README for details"
  }

  command {
    set -euo pipefail

    python3 /tools/create_zarr_optimus.py \
       --cell_metrics ${cell_metrics}\
       --gene_metrics ${gene_metrics}\
       --cell_id ${cell_id}\
       --gene_id  ${gene_id} \
       --output_path_for_zarr "${sample_id}.zarr" \
       --format DirectoryStore \
       --sample_id ${sample_id} \
       --count_matrix ${sparse_count_matrix}

    mkdir zarrout
    # get all the files in the zarr folder in  a list
    a=`find "${sample_id}.zarr"  -type f`
    for f in $a; do
       # replace all / to ! as a work around for now.
       b=`echo $f | tr "/" "\!"`
       mv $f zarrout/$b
    done
  }

  runtime {
    docker: docker
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    preemptible: preemptible
  }

  output {
    Array[File] zarr_output_files = glob("zarrout/*zarr*")
  }
}

