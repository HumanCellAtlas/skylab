import "Optimus.wdl" as target

workflow TestOptimusScientific {
  # Optimus inputs
  Array[File] r1_fastq
  Array[File] r2_fastq

  File whitelist  # 10x genomics cell barcode whitelist for 10x V2
  File tar_star_reference  # star reference
  File annotations_gtf  # gtf containing annotations for gene tagging
  File ref_genome_fasta  # genome fasta file
  String sample_id  # name of sample matching this file, inserted into read group header

  call target.Optimus as target {
    input:
      r1_fastq = r1_fastq,
      r2_fastq = r2_fastq,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
      sample_id = sample_id
  }
  output {
    target.*
  }
}
