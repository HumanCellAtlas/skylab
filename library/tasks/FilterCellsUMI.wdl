task FilterCellsUMI {
    # Input data
    File sparse_count_matrix
    File col_index
    File row_index

    Int n_cells_expected
    String output_file_name = "filter_cells_output.csv"

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-filtercells:0.0.1"
    Int machine_mem_mb = 4400
    Int cpu = 1
    Int disk = 20
    Int preemptible = 3

    meta {
        description: "Detects droplets against empty cells on the basis of the UMI cutoff"
    }

    parameter_meta {
        sparse_count_matrix: "sparse count array in npz format"
        col_index: "sparse count matrix column names in npy format"
        row_index: "sparse count matrix row names in npy format"
        n_cells_expected: "number of cells expected"
        output_file_name: "name of output file"
        cpu: "(optional) the number of cpus to provision for this task"
        disk: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
    }

    command {
        npz2rds.sh -c ${col_index} -r ${row_index} -d ${sparse_count_matrix} -o temp_matrix.rds
        filter_cells.R --input-rds temp_matrix.rds --n_cells_expected ${n_cells_expected} --filter_mode both --output_csv $output_file_name
    }

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        cpu: cpu
        preemptible: preemptible
    }

    output {
        File umi_result = output_file_name
    }
}
