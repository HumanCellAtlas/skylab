workflow SubsetFastqDataset {
    Array[File] bams
    Array[File] fastqs
    String keepregion

    call MergeBams {
        input:
        aligned_bams = bams
    }

    call IndexSortConcatenated {
        input:
        concat_bam = MergeBams.output_concat_bam
    }

    call ExtractReadNames {
        input:
        sorted_bam = IndexSortConcatenated.output_sorted_bam,
        sorted_bai = IndexSortConcatenated.output_sorted_bai,
        subset_region = keepregion
    }

    scatter (fastq in fastqs) {
        call FilterFastq {
            input:
            fastq = fastq,
            keep_read_names_gz = ExtractReadNames.output_kept_reads
        }
    }

    output {
        Array[File] output_filtered_fastq_gz = FilterFastq.output_filtered_fastq_gz
    }
}

task MergeBams {
    Array[File] aligned_bams

    command {
        samtools cat ${sep=' ' aligned_bams} > concat.bam
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
        memory: "8 GB"
        disks: "local-disk 3000 HDD"
        cpu: "2"
    }

    output {
        File output_concat_bam = "concat.bam"
    }
}

task IndexSortConcatenated {
    File concat_bam
    Int samtools_sort_cores = 4
    String samtools_mem_per_core = "5G"

    command {
        samtools sort -@ ${samtools_sort_cores} -m ${samtools_mem_per_core} ${concat_bam} > concat.sorted.bam
        samtools index concat.sorted.bam
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
        memory: "8 GB"
        disks: "local-disk 3000 HDD"
        cpu: "2"
    }

    output {
        File output_sorted_bam = "concat.sorted.bam"
        File output_sorted_bai = "concat.sorted.bam.bai"
    }
}

task ExtractReadNames {
    File sorted_bam
    File sorted_bai
    String subset_region
    Int samtools_sort_cores = 4
    String samtools_mem_per_core = "5G"

    command {
        samtools view -@ ${samtools_sort_cores} -m ${samtools_mem_per_core} ${sorted_bam} ${subset_region} | cut -f 1 | pigz > kept_reads.gz
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
        memory: "28 GB"
        disks: "local-disk 3000 HDD"
        cpu: "8"
    }

    output {
        File output_kept_reads = "kept_reads.gz"
    }
}

task FilterFastq {
    File fastq
    File keep_read_names_gz

    String output_file_name = basename(fastq, ".fastq.gz") + ".filtered.fastq.gz"
    command {
        outfastqgz=${output_file_name}
        filterFastqByReadName.py \
            --in-fastq-gz ${fastq} \
            --out-fastq-gz $outfastqgz \
            --keep-reads-gz ${keep_read_names_gz} \
            --verbose
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
        memory: "28 GB"
        disks: "local-disk 3000 HDD"
        cpu: "8"
    }

    output {
        File output_filtered_fastq_gz = output_file_name
    }
}

task TarFilteredFastq {
    Array[File] filtered_fastqs

    command {
        tar -cvzf output.tar.gz ${sep=' ' filtered_fastqs}
    }

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
        memory: "28 GB"
        disks: "local-disk 3000 HDD"
        cpu: "8"
    }

    output {
        File output_tarred_fastqs = "output.tar.gz"
    }
}