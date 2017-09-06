
# ----------- general questions --------------
# todo are there standard FROM environments we should use? Ubuntu
# todo Array[File] input_fastq <-- how to dereference this?
# sep command in wdl in wdl spec: https://github.com/broadinstitute/wdl/blob/master/SPEC.md
# todo do we have a sense of how we want to publish tools (e.g. a python package containing,
#   among other things, scripts to process data
# --------------------------------------------

task Star {
  File genomic_fastq
  File gtf
  File star_genome  # todo should the gtf just be included in the genome, since they are dependent? ; asked Jishu she said yes

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

  # todo does this execute if local run is requested? set to absurd requirements, see if errors are thrown
  runtime {
    docker: "humancellatlas/star_dev:v1"
    memory: "40 GB"  # todo is input data streamed to the instance? if so, this is too small; jishu answered this
    disks: "local-disk 100 HDD"  # todo is there a distinction between HDD and SSD? yes, but hdd is just as fast at single sequential reads (maybe not paralell reads, ask or benchmark later)
  }
}


task ExtractSubsetReadIndicesInline {
  File bam_file
  Int chromosome

  command <<<
    python <<CODE

    import os
    import json
    import pysam

    # set some constants
    n_aligned = 10000
    n_unaligned = 2000

    with pysam.AlignmentFile(${bam_file}, 'rb') as fin:

    aligned, unaligned = 0, 0  # counters
    chrom_string = str(${chromosome})
    indices = []

    for i, record in enumerate(fin):

        # should be a check that works with all annotation types
        if record.is_unmapped and unaligned < n_unaligned:
            indices.append(i)
            unaligned += 1
        elif not record.is_unmapped and chrom_string in record.rname and aligned < n_aligned:
            indices.append(i)
            aligned += 1

        # check termination condition (we have the requisite number of reads
        if aligned == n_aligned and unaligned == n_unaligned:
            break

    # write indices
    with open(indices.json, 'w') as fout:
        json.dump(indices, fout)

    # warn user if early termination occurred
    if aligned < n_aligned or unaligned < n_unaligned:
        print('Warning: %s: test file construction terminated early. Only %d unaligned '
              'and %d aligned reads were written to %s' %
              (script_name, n_unaligned, n_aligned, indices.json))

    CODE
    >>>

  output {
    File indices_json = "indices.json"
  }

  runtime {
    docker: "ambrosejcarr/python-hca:latest"  # need to make a docker file for this script
    memory: "1 GB"
    disks: "local-disk 100 HDD"
  }

}


task ExtractSubsetReadIndices {
  File bam_file
  Int chromosome

  command {
    extract_testing_reads \
      ${bam_file} \
      indices.json \
      ${chromosome}
  }

  output {
    File output_indices_json = "indices.json"
  }

  runtime {
    docker: "ambrosejcarr/python-hca:latest"
    memory: "1 GB"
    disks: "local-disk 100 HDD"
  }

}


task SubsetFastqFromIndices {
  File indices_json
  File input_fastq_r1
  File input_fastq_r2
  File input_fastq_i1

  # todo is backtick notation the preferred method here?
  command {
    OUTFILES = `subset_fastq_from_indices \
      ${indices_json} \
      ${input_fastq_r1} \
      ${input_fastq_r2} \
      ${input_fastq_i1} \`

    OUTFILES = $(subset_fastq_from_indices \
      ${indices_json} \
      ${input_fastq_r1} \
      ${input_fastq_r2} \
      ${input_fastq_i1})

  }

  runtime {
    docker: "ambrosejcarr/python-hca:latest"
    memory: "2.5 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    Array[String] output_subset_fastqs = read_lines(stdout())
  }

}


task SubsetFastqFromIndicesInline {
  File indices_json
  File input_fastq_r1
  File input_fastq_r2
  File input_fastq_i1

  command <<<
  python <<CODE

  with open(${indices_json}, 'r') as f:
      indices = set(json.load(f))


  # create fastq readers
  fastqs = [${input_fastq_r1}, ${input_fastq_r2}, ${input_fastq_i1}]
  readers = zip(fq.Reader(f) for f in fastqs)

  # create output filenames
  output_filenames = [f.partition('.fastq')[0] + '_subset.fastq.gz' for f in fastqs]

  # open output files
  output_fileobjs = [gzip.open(f, 'w') for f in output_filenames]

  # write records in indices to output files and close files
  try:
      for i, records in enumerate(zip(*readers)):  # bug
          if i in indices:
              for record, fout in zip(records, output_fileobjs):
                  fout.write(record)

  finally:
      for f in output_fileobjs:
          f.close()

  # print to console to enable readlines() of the filenames in output
  for file_ in output_filenames:
      print(file_)

  CODE
  >>>


  runtime {
    docker: "ambrosejcarr/python-hca:latest"
    memory: "2.5 GB"
    disks: "local-disk 1000 HDD"  # todo put this in to see if it breaks cromwell
  }

  output {
    Array[String] output_subset_fastqs = read_lines(stdout())
  }

}


workflow generate_test {

  File input_fastq_r1
  File input_fastq_r2
  File input_fastq_i1
  File gtf
  File star_genome
  Int chromosome

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
      input_fastq_i1 = input_fastq_i1
  }

  output {
    File subset_fastq_r1 = SubsetFastqFromIndices.output_subset_fastqs[0]
    File subset_fastq_r2 = SubsetFastqFromIndices.output_subset_fastqs[1]
    File subset_fastq_i1 = SubsetFastqFromIndices.output_subset_fastqs[2]
  }

}