workflow CellRangerBuildReference {
    File Gencode_GTF
    File Gencode_FASTA
    String Mask_Genes
    String revised_Gencode_GTF = "gencode.v27.primary_assembly.cellranger.annotation.gtf"

    # Runtime attributes
    String memory = "416 GB"
    Int boot_disk_size_gb = 12
    String disk_space = "250"
    Int cpu = 64

    call build_reference {
        input:
        Gencode_GTF = Gencode_GTF,
        Gencode_FASTA = Gencode_FASTA,
        Mask_Genes = Mask_Genes,
        memory = memory,
        boot_disk_size_gb = boot_disk_size_gb,
        disk_space = disk_space,
        cpu = cpu

    }

    output {
        File reference = build_reference.result_reference
    }

}


task build_reference {
    File Gencode_GTF
    File Gencode_FASTA
    String Mask_Genes

    String memory
    Int boot_disk_size_gb
    String disk_space
    Int cpu

    String revised_Gencode_GTF = "gencode.v27.primary_assembly.cellranger.annotation.gtf"
    String reference = "GRCh38"
    String result = "GRCh38_GencodeV27_Primary_CellRanger.tar.gz"

    command {
        set -e

        cellranger mkgtf ${Gencode_GTF} ${revised_Gencode_GTF} ${Mask_Genes}

        cellranger mkref --genome=${reference} --fasta=${Gencode_FASTA} --genes=${revised_Gencode_GTF}  --nthreads=8

        tar cvf ${result} ${reference}
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-cellranger"
        memory: memory
        bootDiskSizeGb: boot_disk_size_gb
        disks: "local-disk " + disk_space + " HDD"
        cpu: cpu
    }

    output {
        File result_reference = "${result}"
    }
}
