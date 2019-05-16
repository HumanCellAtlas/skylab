# test extraction of 10X Genomics v3 chemistry's 12bp UMI barcodes
# (in contrast to 10bp UMI barcodes in v2)
# miniwdl check --path library/tasks --path pipelines/optimus test/optimus/v3umi/TestOptimusV3UMI.wdl

import "Optimus.wdl" as target

workflow TestOptimusV3UMI {
  # Optimus inputs (v3 data with 12bp UMI barcodes)
  Array[File] r1_fastq
  Array[File] r2_fastq
  Array[File] i1_fastq

  File whitelist
  File tar_star_reference
  File annotations_gtf
  File ref_genome_fasta
  String sample_id

  # Run Optimus with tenX_v3_chemistry = true, and check for complete 12bp
  # UMI barcodes in the output BAM file.
  call target.Optimus as targetOptimus12BasePairUMI {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      i1_fastq = i1_fastq,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
      sample_id = sample_id,
      tenX_v3_chemistry = true
  }
  call CheckUMITagLength as check12BasePairUMI {
    input:
      optimus_bam = targetOptimus12BasePairUMI.bam,
      expected_UR_tag_length = 12
  }

  # Run Optimus once more on our v3 data but **without** setting
  # tenX_v3_chemistry = true. Currently, this will "succeed" but produce a
  # garbage output BAM with truncated, 10bp UMI barcodes.
  # In the future, sctools ought to detect & error out in this case, and this
  # test should be adapted accordingly.
  call target.Optimus as targetOptimus10BasePairUMI {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      i1_fastq = i1_fastq,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
      sample_id = sample_id
  }
  call CheckUMITagLength as check10BasePairUMI {
    input:
      optimus_bam = targetOptimus10BasePairUMI.bam,
      expected_UR_tag_length = 10
  }
}

task CheckUMITagLength {
    # verify that all UR tags in bam have the expected value length
    File optimus_bam
    Int expected_UR_tag_length

    # len("UR:Z:") == 5
    Int expected_UR_tag_length_plus_5_base_pairs = expected_UR_tag_length + 5

    command <<<
        set -eux -o pipefail
        L="$(samtools view '${optimus_bam}' \
                | grep -E -o 'UR:Z:[ACGTN]+' \
                | awk '{print length}' \
                | sort -u)"
        declare -r L
        if [ "$L" -ne "${expected_UR_tag_length_plus_5_base_pairs}" ]; then
            >&2 echo "One or more UMI barcode tags had unexpected length other than ~{expected_UR_tag_length} base pairs"
            exit 1
        fi
    >>>

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:1.3.1"
    }
}
