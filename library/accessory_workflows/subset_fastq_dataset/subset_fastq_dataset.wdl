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
      concat_bam = MergeBams.merged_bam
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
        keep_read_names = ExtractReadNames.kept_reads
    }
  }

  output {
    Array[File] filtered_fastqs = FilterFastq.filtered_fastq
  }
}

task MergeBams {
  Array[File] aligned_bams

  Int disk_size = ceil(size(aligned_bams, "GB") * 2) + 10

  command {
    samtools cat ${sep=' ' aligned_bams} > concat.bam
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
    memory: "4 GB"
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
  }

  output {
    File merged_bam = "concat.bam"
  }
}

task IndexSortConcatenated {
  File concat_bam
  Int samtools_sort_cores = 4
  Int samtools_mem_gb_per_core = 5

  Int disk_size = ceil(size(concat_bam, "GB") * 2) + 20

  command {
    samtools sort \
      -@ ${samtools_sort_cores} \
      -m ${samtools_mem_gb_per_core}G \
      ${concat_bam} > concat.sorted.bam

    samtools index concat.sorted.bam
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
    memory: samtools_mem_gb_per_core * samtools_sort_cores
    disks: "local-disk ${disk_size} HDD"
    cpu: ceil((samtools_sort_cores / 2)) + 1
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

  Int disk_size = ceil(size(sorted_bam, "GB")) + 20

  command {
    samtools view \
      ${sorted_bam} \
      ${subset_region} | cut -f 1  > kept_reads
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
    memory: "4 GB"
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
  }

  output {
    File kept_reads = "kept_reads"
  }
}

task FilterFastq {
  File fastq
  File keep_read_names

  String output_file_name = basename(fastq, ".fastq") + ".filtered.fastq"
  Int disk_size = ceil(size(fastq, "GB") * 2) + ceil(size(keep_read_names, "GB")) + 20

  command {
    filterFastqByReadName.py \
      --in-fastq ${fastq} \
      --out-fastq ${output_file_name} \
      --keep-reads ${keep_read_names} \
      --verbose
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
    memory: "2 GB"
    disks: "local-disk ${disk_size} HDD"
    cpu: 1
  }

  output {
    File filtered_fastq = output_file_name
  }
}