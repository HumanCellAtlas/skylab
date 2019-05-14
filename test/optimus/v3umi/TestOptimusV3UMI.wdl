# test extraction of 10X Genomics v3 chemistry's 12bp UMI barcodes
# (in contrast to 10bp UMI barcodes in v2)
# miniwdl check --path library/tasks --path pipelines/optimus test/optimus/v3umi/TestOptimusV3UMI.wdl

import "Optimus.wdl" as target

workflow TestOptimusV3UMI {
  # Optimus inputs (v3 data with 12bp UMI barcodes)
  Array[File] r1_fastq
  Array[File] r2_fastq

  File whitelist
  File tar_star_reference
  File annotations_gtf
  File ref_genome_fasta
  String sample_id

  # First call Optimus **without** setting v3 = true. This should produce
  # a garbage output with incomplete, 10bp UMI barcodes.
  # In the future, sctools ought to detect & error out in this case
  call target.Optimus as target10 {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
      sample_id = sample_id
  }
  call CheckUMITagLength as check10 {
    input:
      optimus_bam = target10.bam,
      expected_UR_tag_length = 10
  }

  # Now try again with v3 = true, and check for complete 12bp UMI barcodes.
  call target.Optimus as target12 {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
      sample_id = sample_id,
      v3 = true
  }
  call CheckUMITagLength as check12 {
    input:
      optimus_bam = target10.bam,
      expected_UR_tag_length = 12
  }
}

task CheckUMITagLength {
    # verify that all UR tags in bam have the expected value length
    File optimus_bam
    Int expected_UR_tag_length
    Int that_plus5 = expected_UR_tag_length + 5 # len("UR:Z:") == 5

    command <<<
        set -eux -o pipefail
        apt-get update && apt-get install -y samtools
        L="$(samtools view '${optimus_bam}' | grep -E -o 'UR:Z:[ACGTN]+' | awk '{print length}' | sort -u)"
        if [ "$L" -ne "${that_plus5}" ]; then
            exit 1
        fi
    >>>

    runtime {
        docker: "ubuntu:18.04"
    }
}
