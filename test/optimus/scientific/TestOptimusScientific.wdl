import "Optimus.wdl" as target
import "ScientificTests.wdl" as tests

workflow TestOptimusScientific {
  Array[Array[File]] fastq_inputs
  File whitelist
  File tar_star_reference
  File annotations_gtf
  File ref_genome_fasta
  String sample_id
  Int expected_number_bam_records

  call target.Optimus as target {
    input:
      fastq_inputs = fastq_inputs,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
  }
      sample_id = sample_id

  call tests.TestBamRecordNumber as test {
    input:
      bam = target.bam,
      expected_records = expected_number_bam_records
  }
}
