task ExtractSubsetReadIndicesInline {
  File bam_file
  Int chromosome

  command <<<
    python <<CODE

    import os
    import sys  # todo debugging
    import json
    import pysam

    # set some constants
    n_aligned = 10000
    n_unaligned = 2000

    with pysam.AlignmentFile('${bam_file}', 'rb') as fin:

        aligned, unaligned = 0, 0  # counters
        chrom_string = str(${chromosome})
        indices = []

        for i, record in enumerate(fin):

            # should be a check that works with all annotation types
            if record.is_unmapped and unaligned < n_unaligned:
                indices.append(i)
                unaligned += 1
            elif not record.is_unmapped and chrom_string in record.reference_name and aligned < n_aligned:
                indices.append(i)
                aligned += 1

            # check termination condition (we have the requisite number of reads
            if aligned == n_aligned and unaligned == n_unaligned:
                break

    # write indices
    with open('indices.json', 'w') as fout:
        json.dump(indices, fout)

    # warn user if early termination occurred
    if aligned < n_aligned or unaligned < n_unaligned:
        print('Warning: %s: test file construction terminated early. Only %d unaligned '
              'and %d aligned reads were written to %s' %
              (script_name, n_unaligned, n_aligned, indices.json))

    CODE
    >>>

  output {
    File output_indices_json = "indices.json"
  }

  runtime {
    docker: "ambrosejcarr/python-hca:latest"  # need to make a docker file for this script
    memory: "1 GB"
    disks: "local-disk 100 HDD"
  }

}


task SubsetFastqFromIndicesInline {
  File indices_json
  File input_fastq_r1
  File input_fastq_r2
  File input_fastq_i1

  command <<<
  python <<CODE

  import json
  import gzip
  import scsequtil.fastq as fq


  # get indices, identify maximum (stopping condition for looping)
  with open('${indices_json}', 'r') as f:
      indices = set(json.load(f))
  max_ind = max(indices)


  # create fastq readers
  fastqs = ['${input_fastq_r1}', '${input_fastq_r2}', '${input_fastq_i1}']
  readers = zip(*[fq.Reader(f) for f in fastqs])

  # create output filenames
  output_filenames = [f.partition('.fastq')[0] + '_subset.fastq.gz' for f in fastqs]

  # open output files
  output_fileobjs = [gzip.open(f, 'wt') for f in output_filenames]

  # write records in indices to output files and close files
  try:
      for i, records in enumerate(readers):
          if i in indices:
              for record, fout in zip(records, output_fileobjs):
                  fout.write(str(record))
          if i > max_ind:
              break

  finally:
      for f in output_fileobjs:
          f.close()

  for f in output_filenames:
      print(f)

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
  File aligned_bam
  File gtf
  File star_genome
  Int chromosome


  call ExtractSubsetReadIndicesInline {
    input:
      bam_file = aligned_bam,
      chromosome = chromosome
  }

  call SubsetFastqFromIndicesInline {
    input:
      indices_json = ExtractSubsetReadIndicesInline.output_indices_json,
      input_fastq_r1 = input_fastq_r1,
      input_fastq_r2 = input_fastq_r2,
      input_fastq_i1 = input_fastq_i1
  }

  output {
    Array[File] output_subset_fastqs = SubsetFastqFromIndicesInline.output_subset_fastqs
  }

}