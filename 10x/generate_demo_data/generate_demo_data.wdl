
task StarAlignSubset {
  File genomic_fastq
  File gtf
  File star_genome
  Int? subset_size

  # note that STAR runThreadN must always equal 1 (default) or the order of the .bam file will be
  # disordered relative to the input data, breaking the expectation of same-ordering necessary
  # for downstream methods
  command {
    # unpack genome reference
    tar -zxvf ${star_genome}

    # truncate the input file
    zcat ${genomic_fastq} | head -n ${default=100000 subset_size} > reads.fastq

    # align reads
    STAR  --readFilesIn reads.fastq \
      --genomeDir ./star \
      --outSAMstrandField intronMotif \
      --genomeLoad NoSharedMemory \
      --sjdbGTFfile ${gtf} \
      --outSAMtype BAM Unsorted  \
      --outSAMunmapped Within \
      --runThreadN 1
  }
  output {
    File output_bam = "Aligned.out.bam"
  }
  runtime {
    docker: "humancellatlas/star_dev:v1"
    memory: "8 GB"  # not used locally
    disks: "local-disk 220 HDD"  # 80G fastq, 32G reference bundle, 80G bam, 15% overflow (30G)
  }
}

task ExtractIndicesSpecificChromosomeAlignments {
  File bam_file
  Int chromosome
  Int? number_alignable_records
  Int? number_unalignable_records

  command <<<
    python3 <<CODE

    import scsequtil.bam as bam
    import json

    sa = bam.SubsetAlignments('${bam_file}')
    alignable, unalignable = sa.indices_by_chromosome(
        ${default=10000 number_alignable_records},
        '${chromosome}',
        ${default=2000 number_unalignable_records})
    with open('indices.json', 'w') as f:
        json.dump(alignable + unalignable, f)

    CODE
  >>>
  output {
    File output_indices_json = "indices.json"
  }
  runtime {
    docker: "ambrosejcarr/python3-scientific:0.2.0"
    memory: "2 GB"
    disks: "local-disk 220 HDD"
  }
}

task SubsetFastqFromIndices {
  File indices_json
  File input_fastq_r1
  File input_fastq_r2
  File input_fastq_i1

  command <<<
    python3 <<CODE

    import json
    import scsequtil.fastq as fq
    import scsequtil.reader as rd
    import gzip

    # get indices
    with open('${indices_json}', 'r') as f:
      indices = set(json.load(f))

    # set fastq inputs
    fastqs = ['${input_fastq_r1}', '${input_fastq_r2}', '${input_fastq_i1}']
    readers = [fq.Reader(f) for f in fastqs]

    # define filenames
    output_filenames = [f.partition('.fastq')[0] + '_subset.fastq.gz' for f in fastqs]

    # open some files
    output_fileobjs = [gzip.open(f, 'wt') for f in output_filenames]

    # write to file
    try:
        for records in rd.zip_readers(*readers, indices=indices):
            for record, fout in zip(records, output_fileobjs):
                fout.write(str(record))
    finally:
        for f in output_fileobjs:
            f.close()

    # write filenames to stdout for cromwell to pick up
    for f in output_filenames:
       print(f)

    CODE
  >>>
  runtime {
    docker: "ambrosejcarr/python3-scientific:0.2.0"
    memory: "2 GB"
    disks: "local-disk 220 HDD"
  }
  output {
    Array[File] output_subset_fastqs = read_lines(stdout())
  }
}

workflow Generate10xDemoData {
  File fastq_r1
  File fastq_r2
  File fastq_i1
  File star_genome
  File gtf
  Int chromosome

  call StarAlignSubset {
    input:
      genomic_fastq = fastq_r2,
      gtf = gtf,
      star_genome = star_genome
  }

  call ExtractIndicesSpecificChromosomeAlignments {
    input:
      bam_file = StarAlignSubset.output_bam,
      chromosome = chromosome
  }

  call SubsetFastqFromIndices {
    input:
      indices_json = ExtractIndicesSpecificChromosomeAlignments.output_indices_json,
      input_fastq_r1 = fastq_r1,
      input_fastq_r2 = fastq_r2,
      input_fastq_i1 = fastq_i1
  }

  output {
    File subset_fastq_r1 = "subset_wf.output_subset_fastqs[0]"
    File subset_fastq_r2 = "subset_wf.output_subset_fastqs[1]"
    File subset_fastq_i1 = "subset_wf.output_subset_fastqs[2]"
  }
}
