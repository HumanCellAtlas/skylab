workflow cellranger {
    String sampleId
    Array[File] fastqs
    String referenceName
    File transcriptomeTarGz
    Int? expectCells
    String diskSpace
    call CellRanger {
        input:
        sampleId = sampleId,
        fastqs = fastqs,
        reference = referenceName,
        transcriptomeTarGz = transcriptomeTarGz,
        expectCells = expectCells,
        diskSpace = diskSpace
   }
}
task CellRanger {
    String sampleId
    Array[File] fastqs
    String reference
    File transcriptomeTarGz
    Int? expectCells
    String diskSpace
    command {
        set -e
        mkdir transcriptome_dir
        tar xf ${transcriptomeTarGz} -C transcriptome_dir --strip-components 1
        ln -s /usr/bin/python3 python
        export PATH=$PATH:.
        python <<CODE
        import os
        from subprocess import call
        dirs = dict()
        for f in ["${sep='","' fastqs}"]:
            dirs.setdefault(os.path.dirname(f), True)
        expect_cells = '${expectCells}'
        call_args = list()
        call_args.append('cellranger')
        call_args.append('count')
        call_args.append('--jobmode=local')
        call_args.append('--transcriptome=transcriptome_dir')
        call_args.append('--id=' + '${sampleId}')
        call_args.append('--fastqs=' + ','.join(list(dirs.keys())))
        call_args.append('--nosecondary')
        call_args.append('--disable-ui')
        if expect_cells is not '':
            call_args.append('--expect-cells=' + str(expect_cells))
        call(call_args)
        CODE
        }
    output {
        File qc = "${sampleId}/outs/metrics_summary.csv"
        File sorted_bam = "${sampleId}/outs/possorted_genome_bam.bam"
        File sorted_bam_index = "${sampleId}/outs/possorted_genome_bam.bam.bai"
        File barcodes = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/barcodes.tsv"
        File genes = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/genes.tsv"
        File matrix = "${sampleId}/outs/filtered_gene_bc_matrices/${reference}/matrix.mtx"
        File filtered_gene_h5 = "${sampleId}/outs/filtered_gene_bc_matrices_h5.h5"
        File raw_gene_h5 = "${sampleId}/outs/raw_gene_bc_matrices_h5.h5"
        File raw_barcodes = "${sampleId}/outs/raw_gene_bc_matrices/${reference}/barcodes.tsv"
        File raw_genes = "${sampleId}/outs/raw_gene_bc_matrices/${reference}/genes.tsv"
        File raw_matrix = "${sampleId}/outs/raw_gene_bc_matrices/${reference}/matrix.mtx"
        File mol_info_h5 = "${sampleId}/outs/molecule_info.h5"
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-cellranger"
        memory: "416 GB"
        bootDiskSizeGb: 12
        disks: "local-disk ${diskSpace} HDD"
        cpu: 64
    }
}

