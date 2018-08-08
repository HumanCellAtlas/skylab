task CellRanger {
    String sample_id
    Array[String]? fastqs
    String? comma_fastqs
    String reference
    File transcriptome_file
    Int? expect_cells
    Boolean? secondary
    Int? force_cells
    Boolean? do_force_cells
    String? chemistry
    Int disk_space
    Int memory
    Int cores
    Int preemptible

    command {
        monitor_script.sh > monitoring.log &

        orchestra_methods.py -c=count \
                            -id=${sample_id} \
                            -cf=${comma_fastqs} \
                            -fs "${sep='" "' fastqs}" \
                            -E=${expect_cells} \
                            -F=${force_cells} \
                            -C=${chemistry} \
                            -S=${secondary} \
                            -tf=${transcriptome_file} \
                            -dfc=${do_force_cells}
        }
    output {
        File barcodes = "results_${sample_id}/outs/filtered_gene_bc_matrices/${reference}/barcodes.tsv"
        File genes = "results_${sample_id}/outs/filtered_gene_bc_matrices/${reference}/genes.tsv"
        File matrix = "results_${sample_id}/outs/filtered_gene_bc_matrices/${reference}/matrix.mtx"
        File qc = "results_${sample_id}/outs/metrics_summary.csv"
        File report = "results_${sample_id}/outs/web_summary.html"
        File sorted_bam = "results_${sample_id}/outs/possorted_genome_bam.bam"
        File sorted_bam_index = "results_${sample_id}/outs/possorted_genome_bam.bam.bai"
        File filtered_gene_h5 = "results_${sample_id}/outs/filtered_gene_bc_matrices_h5.h5"
        File raw_gene_h5 = "results_${sample_id}/outs/raw_gene_bc_matrices_h5.h5"
        File raw_barcodes = "results_${sample_id}/outs/raw_gene_bc_matrices/${reference}/barcodes.tsv"
        File raw_genes = "results_${sample_id}/outs/raw_gene_bc_matrices/${reference}/genes.tsv"
        File raw_matrix = "results_${sample_id}/outs/raw_gene_bc_matrices/${reference}/matrix.mtx"
        File mol_info_h5 = "results_${sample_id}/outs/molecule_info.h5"
        File monitoring_log = "monitoring.log"
    }

 runtime {
         docker: "singlecellportal/scrna-seq_orchestra"
         memory: "${memory} GB"
         bootDiskSizeGb: 12
         disks: "local-disk ${disk_space} HDD"
         cpu: cores
         preemptible: preemptible
     }
}


workflow cellranger {
    meta {
        description: ""
    }

    String sample_id
    Array[String]? fastqs
    String? comma_fastqs
    String reference
    File transcriptome_file
    Boolean secondary = false

    # Cellranger count inputs
    Int? expect_cells
    Int? force_cells
    Boolean? do_force_cells
    String? chemistry

    # Runtime Arguments
    Int? disk_space = 500
    Int? memory = 120
    Int? cores = 32
    Int? preemptible = 2

    parameter_meta {
        sample_id: "Name of sample to run CellRanger count on"
        fastqs: "Array of fastq directories for running decoupled"
        comma_fastqs: "Comma seperated string of fastq directories for running in orchestration"
        reference: "Reference name for count"
        transcriptome_file: "Reference file for count"
        secondary: "Whether to run CellRanger built in secondary analysis"
        expect_cells: "Expected number of recovered cells"
        force_cells: "Number of cells to force the pipeline to use"
        do_force_cells: ""
        chemistry: "Assay configuration"
    }

    call CellRanger {
        input:
        sample_id = sample_id,
        fastqs = fastqs,
        comma_fastqs = comma_fastqs,
        transcriptome_file = transcriptome_file,
        reference = reference,
        secondary = secondary,
        expect_cells = expect_cells,
        disk_space = disk_space,
        force_cells = force_cells,
        do_force_cells = do_force_cells,
        chemistry = chemistry,
        memory = memory,
        cores = cores,
        preemptible = preemptible
    }

    output {
        File barcodes = CellRanger.barcodes
        File genes = CellRanger.genes
        File matrix = CellRanger.matrix
        File qc = CellRanger.qc
        File report = CellRanger.report
        File sorted_bam = CellRanger.sorted_bam
        File sorted_bam_index = CellRanger.sorted_bam_index
        File filtered_gene_h5 = CellRanger.filtered_gene_h5
        File raw_gene_h5 = CellRanger.raw_gene_h5
        File raw_barcodes = CellRanger.raw_barcodes
        File raw_genes = CellRanger.raw_genes
        File raw_matrix = CellRanger.raw_matrix
        File mol_info_h5 = CellRanger.mol_info_h5
        File monitoring_log = CellRanger.monitoring_log
    }
}
