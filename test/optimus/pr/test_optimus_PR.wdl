import "Optimus.wdl" as target
import "ValidateOptimus.wdl" as checker


# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestOptimusPR {

  # output hashes
  String expected_bam_hash
  String expected_matrix_hash
  String expected_gene_metric_hash
  String expected_cell_metric_hash
  String expected_matrix_summary_hash

  # Optimus inputs
  Array[File] r1
  Array[File] r2
  Array[File]? i1

  File whitelist  # 10x genomics cell barcode whitelist for 10x V2
  File tar_star_reference  # star reference
  File annotations_gtf  # gtf containing annotations for gene tagging
  File ref_genome_fasta  # genome fasta file
  String sample_id  # name of sample matching this file, inserted into read group header

  call target.Optimus as target {
    input:
      r1 = r1,
      r2 = r2,
      i1 = i1,
      whitelist = whitelist,
      tar_star_reference = tar_star_reference,
      annotations_gtf = annotations_gtf,
      ref_genome_fasta = ref_genome_fasta,
      sample_id = sample_id
  }

  call checker.ValidateOptimus as checker {
    input:
      bam = target.bam,
      matrix = target.matrix,
      gene_metrics = target.gene_metrics,
      cell_metrics = target.cell_metrics,
      matrix_summary = target.matrix_summary,
      expected_bam_hash = expected_bam_hash,
      expected_matrix_hash = expected_matrix_hash,
      expected_cell_metric_hash = expected_cell_metric_hash,
      expected_gene_metric_hash = expected_gene_metric_hash,
      expected_matrix_summary_hash = expected_matrix_summary_hash
  }

}
