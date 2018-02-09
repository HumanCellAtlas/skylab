
task StarAlignBamSingleEnd {
  File bam_input  # unaligned bam file containing genomic sequence, tagged with barcode information
  File tar_star_reference  # star reference tarball

  # estimate input requirements for star using UBams:
  # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
  Int input_size = ceil(size(bam_input, "G"))
  Int reference_size = ceil(size(tar_star_reference, "G"))
  Int estimated_disk_required =  ceil(input_size * 2.2 + reference_size * 2)
  Int estimated_memory_required = reference_size + 6

  command {
    # prepare reference
    mkdir genome_reference
    tar -xf "${tar_star_reference}" -C genome_reference --strip-components 1
    rm "${tar_star_reference}"

    STAR \
      --runMode alignReads \
      --runThreadN $(nproc) \
      --genomeDir genome_reference \
      --readFilesIn "${bam_input}" \
      --outSAMtype BAM Unsorted \
      --outSAMattributes All \
      --outFilterMultimapNmax 1 \
      --outSAMunmapped Within \
      --outSAMprimaryFlag AllBestScore \
      --readFilesType SAM SE \
      --readFilesCommand samtools view -h \
      --runRNGseed 777
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-40ead6e"
    cpu: 16
    memory: "${estimated_memory_required} GB"
    disks: "local-disk ${estimated_disk_required} SSD"
  }

  output {
    File bam_output = "Aligned.out.bam"
    File alignment_log = "Log.final.out"
  }

}
