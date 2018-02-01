
task StarAlignFastqSingleEnd {
  File fastq_input  # fastq file containing genomic sequence
  File tar_star_reference  # gzipped star reference tarball

  # estimate input requirements for star using uBams:
  # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
  Int input_size = ceil(size(fastq_input, "G"))
  Int reference_size = ceil(size(tar_star_reference, "G"))
  Int estimated_disk_required =  ceil(input_size * 2.2 + reference_size * 2)

  command {
    # prepare reference
    mkdir genome_reference && tar -xf "${tar_star_reference}" -C genome_reference --strip-components 1
    rm "${tar_star_reference}"

    STAR \
      --runMode alignReads \
      --runThreadN $(nproc) \
      --genomeDir genome_reference \
      --readFilesIn "${fastq_input}" \
      --outSAMtype BAM Unsorted \
      --readFilesCommand zcat

  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-40ead6e"
    cpu: 16
    memory: "30 GB"
    disks: "local-disk ${estimated_disk_required} SSD"
  }

  output {
    File bam = "Aligned.out.bam"
    File alignment_log = "Log.final.out"
  }

}
