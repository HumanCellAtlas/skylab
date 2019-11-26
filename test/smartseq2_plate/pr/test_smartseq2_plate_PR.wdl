import "SmartSeq2Plate.wdl" as target_wdl
import "ValidateSmartSeq2Plate.wdl" as checker_wdl

# this task will be run by the jenkins script that gets executed on our PRs.
workflow TestSmartSeq2SingleCellPR {

  # expected hashes of target_workflow outputs
  String expected_core_QC_hash
  String expected_qc_tabls_hash
  String expected_gene_matrix_hash
  String expected_isoform_matrix_hash

  # SS2 inputs
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat
  String hisat2_ref_name
  String hisat2_ref_trans_name
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index
  String stranded

  # Plate information and input files
  String file_prefix
  Array[String] input_file_names
  String batch_id

  call target_wdl.RunSmartSeq2ByPlate as target_workflow {
    input:
      genome_ref_fasta = genome_ref_fasta,
      rrna_intervals = rrna_intervals,
      gene_ref_flat = gene_ref_flat,
      hisat2_ref_name = hisat2_ref_name,
      hisat2_ref_trans_name = hisat2_ref_trans_name,
      hisat2_ref_index = hisat2_ref_index,
      hisat2_ref_trans_index = hisat2_ref_trans_index,
      rsem_ref_index = rsem_ref_index,
      stranded = stranded,
      file_prefix = file_prefix,
      input_file_names = input_file_names,
      batch_id = batch_id
  }

  call checker_wdl.ValidateSmartSeq2Plate as checker_workflow {
    input:
     core_QC = target_workflow.core_QC,
     qc_tabls = target_workflow.qc_tabls,
     gene_matrix = target_workflow.gene_matrix,
     isoform_matrix = target_workflow.isoform_matrix,
     expected_core_QC_hash = expected_core_QC_hash,
     expected_qc_tabls_hash = expected_qc_tabls_hash,
     expected_gene_matrix_hash = expected_gene_matrix_hash,
     expected_isoform_matrix_hash =expected_isoform_matrix_hash
  }

}
