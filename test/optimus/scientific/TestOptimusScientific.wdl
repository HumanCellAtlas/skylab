import "Optimus.wdl" as target

workflow TestOptimusScientific {

  meta {
    description: "This workflow runs the Optimus pipeline on synthetic testing data. #TODO it should also test this data against the output of the Cell Ranger reference pipeline"
  }

  # Optimus inputs
  Array[File] r1_fastq
  Array[File] r2_fastq

  File whitelist
  File tar_star_reference
  File annotations_gtf
  File ref_genome_fasta
  String sample_id

  parameter_meta {
    r1_fastq: "The forward (barcode containing) read of a 10x genomics RNA v2 sequencing experiment"
    r2_fastq: "The reverse (genomic sequence containing) read of a 10x genomics RNA v2 sequencing experiment"
    whitelist: "A list of expected cell barcodes"
    annotations_gtf: "A gtf containing annotations for gene tagging"
    ref_genome_fastq: "The genome reference"
    sample_id: "Name of the sample"
  }

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
