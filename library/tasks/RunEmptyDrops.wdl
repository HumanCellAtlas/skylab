task RunEmptyDrops {
    # Input data
    File sparse_count_matrix
    File col_index
    File row_index

    # emptyDrops Params
    Float niters = 10000
    Float fdr_cutoff = 0.10

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-dropletutils:0.1.0"
    Int machine_mem_mb = 4000
    Int cpu = 1
    Int disk = 20
    Int preemptible = 0

    meta {
        description: "Runs empty drops on the count matrix and calls cells on the  basis of the provided cutoff value"
    }

    parameter_meta {
        sparse_count_matrix: "sparse count array in npz format"
        col_index: "sparse count matrix column names in npy format"
        row_index: "sparse count matrix row names in npy format"
        cpu: "(optional) the number of cpus to provision for this task"
        disk: "(optional) the amount of disk space (GB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command {
        npz2rds.sh -c ${col_index} -r ${row_index} -d ${sparse_count_matrix} -o temp_matrix.rds
        emptyDropsWrapper.R --transpose --verbose --input-rds temp_matrix.rds --output-csv empty_drops_result.csv --fdr-cutoff ${fdr_cutoff} --emptydrops-niters ${niters}
    }

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File empty_drops_result = "empty_drops_result.csv"
    }
}
