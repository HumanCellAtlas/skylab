workflow SubsetFastqDataset {
    String alignment_location
    String fastq_location
    String keepregion

    call DownloadMergeBams {
        input:
        star_align_bucket = alignment_location
    }

    call IndexSortConcatenated {
        input:
        concat_bam = DownloadMergeBams.output_concat_bam
    }

    call ExtractReadNames {
        input:
        sorted_bam = IndexSortConcatenated.output_sorted_bam,
        sorted_bai = IndexSortConcatenated.output_sorted_bai,
        subset_region = keepregion
    }

    call DownloadFastqs {
        input:
        in_fastq_bucket = fastq_location
    }

    scatter (fastq in DownloadFastqs.output_fastqs) {
        call FilterFastq {
            input:
            fastq = fastq,
            keep_read_names_gz = ExtractReadNames.output_kept_reads
        }
    }

    call TarFilteredFastq {
        input:
        filtered_fastqs = FilterFastq.output_filtered_fastq_gz
    }

    output {
        File output_tar = TarFilteredFastq.output_tarred_fastqs
    }
}

task DownloadMergeBams {
    String star_align_bucket

    command {
        mkdir starAlignOutputs
        gsutil -m cp -r ${star_align_bucket} starAlignOutputs/
        samtools cat `find ./starAlignOutputs/ -name '*.bam' -printf '%p '` > concat.bam
    }

     runtime {
       docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
       memory: "28 GB"
       disks: "local-disk 3000 HDD"
       cpu: "8"
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
       memory: "28 GB"
       disks: "local-disk 3000 HDD"
       cpu: "8"
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

    command {
        samtools view ${sorted_bam} ${subset_region} | cut -f 1 | pigz > kept_reads.gz
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

task DownloadFastqs {
    String in_fastq_bucket

    command {
        gsutil -m cp -r ${in_fastq_bucket} ./
    }

     runtime {
       docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
       memory: "28 GB"
       disks: "local-disk 3000 HDD"
       cpu: "8"
     }

    output {
        Array[File] output_fastqs = glob("*.fastq.gz")
    }
}

task FilterFastq {
    File fastq
    File keep_read_names_gz

    command {
        outfastqgz=`basename ${fastq}`
        filterFastqBtReadName.py \
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
        File output_filtered_fastq_gz = "output.gz"
    }
}

task TarFilteredFastq {
    Array[File] filtered_fastqs

    command {
        tar cvzf output.tar.gz ${filtered_fastqs}
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