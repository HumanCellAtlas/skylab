workflow CellRanger {
    meta {
        description: "Analyze 3' single-cell RNA-seq data using the 10X Genomics Cellranger pipeline."
    }
    # version of this pipeline
    String version = "cellranger_v1.0.2"

    String sample_id
    Array[File] fastqs
    String reference_name
    File transcriptome_tar_gz
    Int? expect_cells

    # Runtime attributes
    String memory = "416 GB"
    Int boot_disk_size_gb = 12
    String disk_space = "400"
    Int cpu = 64
    Int max_retries = 0

    parameter_meta {
        sample_id: "Name of sample to run CellRanger count on"
        fastqs: "Array of fastq files"
        reference_name: "Reference name for count"
        transcriptome_tar_gz: "CellRanger-compatible transcriptome reference (can be generated with cellranger mkref)"
        expect_cells: "Expected number of recovered cells (defaults to 3000)"
        memory: "The minimum amount of RAM to use for the Cromwell VM"
        boot_disk_size_gb: "Size of disk (GB) where the docker image is booted by the Cromwell VM"
        disk_space: "Amount of disk space (GB) to allocate to the Cromwell VM"
        cpu: "The minimum number of cores to use for the Cromwell VM"
        max_retries: "(optional) retry this number of times if task fails -- use with caution, see skylab README for details"
    }

    call cellranger_count {
        input:
        sample_id = sample_id,
        fastqs = fastqs,
        reference = reference_name,
        transcriptome_tar_gz = transcriptome_tar_gz,
        expect_cells = expect_cells,
        memory = memory,
        boot_disk_size_gb = boot_disk_size_gb,
        disk_space = disk_space,
        cpu = cpu,
        max_retries = max_retries
   }

   output {
       # version of this pipeline
       String pipeline_version = version

       File qc = cellranger_count.qc
       File sorted_bam = cellranger_count.sorted_bam
       File sorted_bam_index = cellranger_count.sorted_bam_index
       File barcodes = cellranger_count.barcodes
       File genes = cellranger_count.genes
       File matrix = cellranger_count.matrix
       File filtered_gene_h5 = cellranger_count.filtered_gene_h5
       File raw_gene_h5 = cellranger_count.raw_gene_h5
       File raw_barcodes = cellranger_count.raw_barcodes
       File raw_genes = cellranger_count.raw_genes
       File raw_matrix = cellranger_count.raw_matrix
       File mol_info_h5 = cellranger_count.mol_info_h5
       File web_summary = cellranger_count.web_summary
   }
}

task cellranger_count {
    String sample_id
    Array[File] fastqs
    String reference
    File transcriptome_tar_gz
    Int? expect_cells
    String memory
    Int boot_disk_size_gb
    String disk_space
    Int cpu
    Int max_retries = 0

    command {
        set -e

        # The cellranger count transcriptome parameter requires a directory named transcriptome_dir
        # that contains fasta/genome.fa and genes/genes.gtf
        mkdir transcriptome_dir
        tar xf ${transcriptome_tar_gz} -C transcriptome_dir --strip-components 1

        python <<CODE
        import os
        from subprocess import check_call

        # Get paths to fastq directories
        fastq_dirs = set([os.path.dirname(f) for f in "${sep='","' fastqs}"])

        # Call cellranger count
        call_args = list()
        call_args.append('cellranger')
        call_args.append('count')
        call_args.append('--jobmode=local')
        call_args.append('--transcriptome=transcriptome_dir')
        call_args.append('--id=' + '${sample_id}')
        call_args.append('--fastqs=' + ','.join(list(fastq_dirs)))
        call_args.append('--nosecondary')
        call_args.append('--disable-ui')

        expect_cells = '${expect_cells}'
        if expect_cells:
            call_args.append('--expect-cells=' + str(expect_cells))
        check_call(call_args)

        CODE
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-cellranger:v1.0.0"
        memory: memory
        bootDiskSizeGb: boot_disk_size_gb
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
        max_retries: max_retries
    }

    output {
        File qc = "${sample_id}/outs/metrics_summary.csv"
        File sorted_bam = "${sample_id}/outs/possorted_genome_bam.bam"
        File sorted_bam_index = "${sample_id}/outs/possorted_genome_bam.bam.bai"
        File barcodes = "${sample_id}/outs/filtered_gene_bc_matrices/${reference}/barcodes.tsv"
        File genes = "${sample_id}/outs/filtered_gene_bc_matrices/${reference}/genes.tsv"
        File matrix = "${sample_id}/outs/filtered_gene_bc_matrices/${reference}/matrix.mtx"
        File filtered_gene_h5 = "${sample_id}/outs/filtered_gene_bc_matrices_h5.h5"
        File raw_gene_h5 = "${sample_id}/outs/raw_gene_bc_matrices_h5.h5"
        File raw_barcodes = "${sample_id}/outs/raw_gene_bc_matrices/${reference}/barcodes.tsv"
        File raw_genes = "${sample_id}/outs/raw_gene_bc_matrices/${reference}/genes.tsv"
        File raw_matrix = "${sample_id}/outs/raw_gene_bc_matrices/${reference}/matrix.mtx"
        File mol_info_h5 = "${sample_id}/outs/molecule_info.h5"
        File web_summary = "${sample_id}/outs/web_summary.html"
    }
}
