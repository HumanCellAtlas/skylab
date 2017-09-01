import "prepare_and_analyze_ss2_single_sample.wdl" as prep_and_analyze
import "submit_ss2_single_sample.wdl" as submit

workflow WrapperSs2RsemSingleSample {
  String bundle_uuid
  String bundle_version

  File gtf
  File ref_fasta
  File rrna_interval
  File ref_flat
  String star_genome
  String rsem_genome

  call prep_and_analyze.PrepareAndAnalyzeSs2RsemSingleSample as prep_analyze {
    input:
      bundle_uuid = bundle_uuid,
      bundle_version = bundle_version,
      gtf = gtf,
      ref_fasta = ref_fasta,
      rrna_interval = rrna_interval,
      ref_flat = ref_flat,
      star_genome = star_genome,
      rsem_genome = rsem_genome
  }

  call submit.SubmitSs2RsemSingleSample {
    input:
      bam_file = prep_analyze.bam_file,
      bam_trans = prep_analyze.bam_trans,
      rna_metrics = prep_analyze.rna_metrics,
      aln_metrics = prep_analyze.aln_metrics,
      rsem_gene_results = prep_analyze.rsem_gene_results,
      rsem_isoform_results = prep_analyze.rsem_isoform_results,
      rsem_gene_count = prep_analyze.rsem_gene_count,
      gene_unique_counts = prep_analyze.gene_unique_counts,
      exon_unique_counts = prep_analyze.exon_unique_counts,
      transcript_unique_counts = prep_analyze.transcript_unique_counts
  }
}
