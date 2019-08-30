task SmartSeq2ZarrConversion {

  #runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.10"
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
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
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
    memory: "18 GiB"
    disks: "local-disk 100 HDD"
    preemptible: preemptible
  }

  output {
    Array[File] zarr_output_files = glob("zarrout/*zarr*")
  }
}


task OptimusZarrConversion {
  #runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-zarr-output:0.0.1"

  # name of the sample
  String sample_id
  # gene annotation file in GTF format
  File annotation_file
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
    machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"
    cpu: "(optional) the number of cpus to provision for this task"
    disk: "(optional) the amount of disk space (GiB) to provision for this task"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  command {
    set -euo pipefail

    python3 /tools/create_zarr_optimus.py \
       --empty_drops_file ${empty_drops_result} \
       --annotation_file ${annotation_file}\
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
    memory: "18 GiB"
    disks: "local-disk 100 HDD"
    preemptible: preemptible
  }

  output {
    Array[File] zarr_output_files = glob("zarrout/*zarr*")
  }
}

task OptimusZarrToLoom {
     # runtime values
     String docker = "quay.io/humancellatlas/zarr-to-loom:0.0.1"

     Array[File] zarr_files

     Int preemptible = 3
     Int cpu = 1

     meta {
         description: "This task converts the Optimus Zarr output into a loom file"
     }

     parameter_meta {
         machine_mem_mb: "(optional) the amount of memory in (MiB) to provision for this task"
	 cpu: "(optional) the number of cpus to provision for this task"
	 preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non-preemptible machine"
    }
    
    command {
       set -euo pipefail

       mkdir packed_zarr
       mv $zarr_files packed_zarr
       mkdir unpacked_zarr
       ./unpackZARR.sh -i packed_zarr -o unpacked_zarr
       ./optimus_zarr_to_loom.py --input-zarr unpacked_zarr --output-loom output.loom
    }

    runtime {nn
        docker: docker
	cpu: 1
	memory: "10 GiB"
	disks: "local-disk 100 HDD"
	preemptible: preemptible
    }

    output {
       File loom_output = "output.loom"
    }
}