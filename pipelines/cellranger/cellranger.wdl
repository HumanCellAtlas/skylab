workflow CellRanger {
    meta {
        description: "Analyze 3' single-cell RNA-seq data using the 10X Genomics Cellranger pipeline."
    }

    String sample_id
    Array[File] fastqs
    String reference_name
    File transcriptome_tar_gz
    Int? expect_cells

    # Runtime attributes
    String memory = "416 GB"
    Int boot_disk_size_gb = 12
    String disk_space = "250"
    Int cpu = 64

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
        cpu = cpu
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

    command {
        set -e

        # The cellranger count transcriptome parameter requires a directory named transcriptome_dir
        # that contains fasta/genome.fa and genes/genes.gtf
        mkdir transcriptome_dir
        tar xf ${transcriptome_tar_gz} -C transcriptome_dir --strip-components 1
        ln -s /usr/bin/python3 python
        export PATH=$PATH:.

        python <<CODE
        import os
        from subprocess import call

        # Get paths to fastq directories
        dirs = dict()
        for f in ["${sep='","' fastqs}"]:
            dirs.setdefault(os.path.dirname(f), True)

        # Call cellranger count
        call_args = list()
        call_args.append('cellranger')
        call_args.append('count')
        call_args.append('--jobmode=local')
        call_args.append('--transcriptome=transcriptome_dir')
        call_args.append('--id=' + '${sample_id}')
        call_args.append('--fastqs=' + ','.join(list(dirs.keys())))
        call_args.append('--nosecondary')
        call_args.append('--disable-ui')

        expect_cells = '${expect_cells}'
        if expect_cells:
            call_args.append('--expect-cells=' + str(expect_cells))
        call(call_args)

        CODE
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-cellranger"
        memory: memory
        bootDiskSizeGb: boot_disk_size_gb
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
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
    }
}
