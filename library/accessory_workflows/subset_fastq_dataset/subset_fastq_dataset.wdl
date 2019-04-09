version 1.0

workflow SubsetFastqDataset {

  input {
    Array[File] bams
    Array[File] fastqs
    String keepregion
  }

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

  input {
    Array[File] aligned_bams
  }

  Int disk_size = ceil(size(aligned_bams, "GiB") * 10) + 100

  command {
    samtools cat ~{sep=' ' aligned_bams} > concat.bam
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
    memory: "4 GiB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
  }

  output {
    File merged_bam = "concat.bam"
  }
}

task IndexSortConcatenated {

  input {
    File concat_bam
    Int samtools_sort_cores = 4
    Int samtools_mem_gb_per_core = 5
  }

  Int disk_size = ceil(size(concat_bam, "GiB") * 2) + 20

  command {
    samtools sort \
      -@ ~{samtools_sort_cores} \
      -m ~{samtools_mem_gb_per_core}G \
      ~{concat_bam} > concat.sorted.bam

    samtools index concat.sorted.bam
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
    memory: "~{samtools_mem_gb_per_core * samtools_sort_cores} GiB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: ceil((samtools_sort_cores / 2)) + 1
  }

  output {
    File output_sorted_bam = "concat.sorted.bam"
    File output_sorted_bai = "concat.sorted.bam.bai"
  }
}

task ExtractReadNames {

  input {
    File sorted_bam
    File sorted_bai
    String subset_region
  }

  Int disk_size = ceil(size(sorted_bam, "GiB")) + 20

  command {
    samtools view \
      ~{sorted_bam} \
      ~{subset_region} | cut -f 1  > kept_reads
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
    memory: "4 GiB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
  }

  output {
    File kept_reads = "kept_reads"
  }
}

task FilterFastq {

  input {
    File fastq
    File keep_read_names
  }

  String output_file_name = basename(fastq, ".fastq.gz") + ".filtered.fastq.gz"
  Int disk_size = ceil(size(fastq, "GiB") * 2) + ceil(size(keep_read_names, "GiB")) + 20

  command {
    filterFastqByReadName.py \
      --in-fastq-gz ~{fastq} \
      --out-fastq-gz ~{output_file_name} \
      --keep-reads ~{keep_read_names} \
      --verbose
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
    memory: "2 GiB"
    disks: "local-disk ~{disk_size} HDD"
    cpu: 1
  }

  output {
    File filtered_fastq = output_file_name
  }
}