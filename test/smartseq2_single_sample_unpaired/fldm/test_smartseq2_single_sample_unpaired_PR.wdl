import "SmartSeq2SingleSampleUnpaired.wdl" as target_wdl
import "ValidateSmartSeq2SingleCellUnpairedFldm.wdl" as checker_wdl

# this task will be run by the jenkins script that gets executed on our PRs.
workflow TestSmartSeq2SingleCellFluidigmUnpairedPR {

  # expected hashes of target_workflow outputs
  String expected_gene_counts_hash

  # SS2 inputs
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index
  String hisat2_ref_name
  String hisat2_ref_trans_name
  String stranded
  String sample_name
  String output_name
  File fastq

  call target_wdl.SmartSeq2SingleCellUnpaired as target_workflow {
    input:
      genome_ref_fasta = genome_ref_fasta,
      rrna_intervals = rrna_intervals,
      gene_ref_flat = gene_ref_flat,
      hisat2_ref_index = hisat2_ref_index,
      hisat2_ref_trans_index = hisat2_ref_trans_index,
      rsem_ref_index = rsem_ref_index,
      hisat2_ref_name = hisat2_ref_name,
      hisat2_ref_trans_name = hisat2_ref_trans_name,
      stranded = stranded,
      sample_name = sample_name,
      output_name = output_name,
      fastq = fastq
  }

  call checker_wdl.ValidateSmartSeq2SingleCellUnpairedFldm as checker_workflow {
    input:
     alignment_summary_metrics = target_workflow.alignment_summary_metrics,
     gene_counts = target_workflow.rsem_gene_results,
     expected_gene_counts_hash = expected_gene_counts_hash
  }

}
