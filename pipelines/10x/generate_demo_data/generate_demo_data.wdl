
task StarAlignSubset {
  File genomic_fastq
  File gtf
  File star_genome
  Int? subset_size = 1000000

  # note that STAR runThreadN must always equal 1 (default) or the order of the .bam file will be
  # disordered relative to the input data, breaking the expectation of same-ordering necessary
  # for downstream methods
  command {
    # unpack genome reference
    tar -zxvf ${star_genome}

    # truncate the input file
    zcat ${genomic_fastq} | head -n ${subset_size} > reads.fastq

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
    docker: "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-1.0.0"
    memory: "8 GB"  # not used locally
    disks: "local-disk 220 HDD"  # 80G fastq, 32G reference bundle, 80G bam, 15% overflow (30G)
  }
}

task ExtractIndicesSpecificChromosomeAlignments {
  File bam_file
  Int chromosome
  Int? number_alignable_records = 10000
  Int? number_unalignable_records = 2000

  command <<<
    python3 <<CODE

    import scsequtil.bam as bam
    import json

    sa = bam.SubsetAlignments('${bam_file}')
    alignable, unalignable = sa.indices_by_chromosome(
        ${number_alignable_records},
        '${chromosome}',
        ${number_unalignable_records})
    with open('indices.json', 'w') as f:
        json.dump(alignable + unalignable, f)

    CODE
  >>>
  output {
    File output_indices_json = "indices.json"
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.5"
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

    import os
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
    filenames_nopath = [os.path.split(f)[1] for f in fastqs]
    output_filenames = [f.partition('.fastq')[0] + '_subset.fastq.gz' for f in filenames_nopath]

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

    CODE
  >>>
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.5"
    memory: "2 GB"
    disks: "local-disk 220 HDD"
  }
  output {
    Array[File] output_subset_fastqs = glob("./*_subset.fastq.gz")
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
    File subset_fastq_r1 = SubsetFastqFromIndices.output_subset_fastqs[0]
    File subset_fastq_r2 = SubsetFastqFromIndices.output_subset_fastqs[1]
    File subset_fastq_i1 = SubsetFastqFromIndices.output_subset_fastqs[2]
  }
}
