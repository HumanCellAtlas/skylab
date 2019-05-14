import "bulk_rna_pipeline.wdl" as target
import "ValidateBulkRna.wdl" as checker


# this workflow will be run by the jenkins script that gets executed by PRs.
workflow TestBulkRnaPR {

  # output hashes
  String expected_anno_flagstat_hash
  String expected_annobam_hash
  String expected_genome_flagstat_hash
  String expected_genomebam_hash
  String expected_quants_hash
  String expected_mad_qc_metrics_hash
  String expected_rna_qc_metrics_hash
  String expected_rsem_genes_detected_hash

  # Optimus inputs
  Array[File] fastqs_R1
  Array[File] fastqs_R2

  File align_index
  File rsem_index
  File kallisto_index
  File chrom_sizes
  File rna_qc_tr_id_to_gene_type_tsv

  call target.rna as target {
    input:
      endedness = "paired",
      fastqs_R1 = fastqs_R1,
      fastqs_R2 = fastqs_R2,
      aligner = "star",
      align_index = align_index,
      rsem_index = rsem_index,
      kallisto_index = kallisto_index,
      bamroot = "PE_stranded",
      strandedness = "stranded",
      strandedness_direction = "reverse",
      chrom_sizes = chrom_sizes,
      align_ncpus = 2,
      align_ramGB = 4,
      rsem_ncpus = 2,
      rsem_ramGB = 4,
      kallisto_number_of_threads = 2,
      kallisto_ramGB = 4,
      rna_qc_tr_id_to_gene_type_tsv = rna_qc_tr_id_to_gene_type_tsv,
      bam_to_signals_ncpus = 1,
      bam_to_signals_ramGB = 2,
      align_disk = "local-disk 20 HDD",
      kallisto_disk = "local-disk 20 HDD",
      rna_qc_disk = "local-disk 20 HDD",
      bam_to_signals_disk = "local-disk 20 HDD",
      mad_qc_disk = "local-disk 20 HDD",
      rsem_disk = "local-disk 20 HDD"
  }

  call checker.ValidateBulkRna as checker {
    input:
      anno_flagstat = target.anno_flagstat[0],
      annobam = target.annobam[0],
      genome_flagstat = target[0],
      genomebam = target[0],
      bam_to_signals_all = target.bam_to_signals_all[0],
      quants = target.quants[0],
      mad_qc_metrics = target.mad_qc_metrics,
      rna_qc_metrics = target.rna_qc_metrics[0],
      rsem_genes_detected = target.rsem_genes_detected[0],
      expected_anno_flagstat_hash = expected_anno_flagstat_hash,
      expected_annobam_hash = expected_annobam_hash,
      expected_genome_flagstat_hash = expected_genome_flagstat_hash,
      expected_genomebam_hash = expected_genomebam_hash,
      expected_quants_hash = expected_quants_hash,
      expected_mad_qc_metrics_hash = expected_mad_qc_metrics_hash,
      expected_rna_qc_metrics_hash = expected_rna_qc_metrics_hash,
      expected_rsem_genes_detected_hash = expected_rsem_genes_detected_hash
  }

}
