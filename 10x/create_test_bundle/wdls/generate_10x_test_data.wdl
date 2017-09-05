
# ----------- general questions --------------
# todo are there standard FROM environments we should use?
# todo Array[File] input_fastq <-- how to dereference this?
# todo do we have a sense of how we want to publish tools (e.g. a python package containing,
#   among other things, scripts to process data
# --------------------------------------------

task Star {
  File genomic_fastq
  File gtf
  File star_genome  # todo should the gtf just be included in the genome, since they are dependent?

  # note that if runThreadN must always equal 1 or the order will be modified
  # could pipe zcat to head in readFilesCommand for speed? need only align about ~ 10m reads
  command {
    tar -zxvf ${star_genome}
    STAR  --readFilesIn ${genomic_fastq} \
      --genomeDir ./star \
      --quantMode TranscriptomeSAM \
      --outSAMstrandField intronMotif \
      --genomeLoad NoSharedMemory \
      --sjdbGTFfile ${gtf} \
      --readFilesCommand "zcat" \
      --outSAMtype BAM Unsorted  \
      --outSAMunmapped Within \
      --limitBAMsortRAM 30000000000
  }
  output {
    File output_bam = "Aligned.out.bam"
  }

  # todo does this execute if local run is requested?
  runtime {
    docker: "humancellatlas/star_dev:v1"
    memory: "40 GB"  # todo is input data streamed to the instance? if so, this is too small
    disks: "local-disk 100 HDD"  # todo is there a distinction between HDD and SSD?
  }
}


# todo do I need to make a repository for all these python scripts/tools?
task ExtractSubsetReadIndices {
  File bam_file
  File indices_json
  String output_fastq_prefix
  Int chromosome

  # todo
  command {
    extract_testing_reads \
      ${bam_file} \
      ${output_fastq_prefix}_indices.json \
      ${chromosome}
  }

  output {
    File output_indices_json = "${output_fastq_prefix}_indices.json"
  }

  runtime {
    docker: ""  # need to make a docker file for this script
    memory: "1 GB"
    disks: "local-disk 100 HDD"
  }

}

# todo if passing a filename for a future file, do I pass as file or string?
task SubsetFastqFromIndices {
  File indices_json
  File input_fastq_r1
  File input_fastq_r2
  File input_fastq_i1
  String output_fastq_prefix

  # todo is backtick notation the preferred method here?
  command {
    OUTFILES = `subset_fastq_from_indices \
      ${indices_json} \
      ${input_fastq_r1} \
      ${input_fastq_r2} \
      ${input_fastq_i1} \`

  }

  runtime {
    docker: ""  # ??
    memory: "2.5 GB"  # todo if I give picard 2gb, is this enough overflow ram?
    disks: "local-disk 100 HDD"
  }

  output {
    File subset_fastq_r1 = "${OUTFILES[0]}"
    File subset_fastq_r2 = "${OUTFILES[1]}"
    File subset_fastq_i1 = "${OUTFILES[2]}"
  }

}

workflow generate_test {

  File input_fastq_r1
  File input_fastq_r2
  File input_fastq_i1
  File gtf
  File star_genome
  Int chromosome
  String output_fastq_prefix

  call Star {
    input:
      genomic_fastq = input_fastq_r2,
      gtf = gtf,
      star_genome = star_genome
  }

  call ExtractSubsetReadIndices {
    input:
      bam_file = Star.output_bam,
      chromosome = chromosome
  }

  call SubsetFastqFromIndices {
    input:
      indices_json = ExtractSubsetReadIndices.output_indices_json,
      input_fastq_r1 = input_fastq_r1,
      input_fastq_r2 = input_fastq_r2,
      input_fastq_i1 = input_fastq_i1,
      output_fastq_prefix = output_fastq_prefix
  }

  output {
    File subset_fastq_r1 = SubsetFastqFromIndices.output_subset_fastq_r1
    File subset_fastq_r2 = SubsetFastqFromIndices.output_subset_fastq_r2
    File subset_fastq_i1 = SubsetFastqFromIndices.output_subset_fastq_i1
  }

}