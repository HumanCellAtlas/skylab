import "HISAT2QualityControlPipeline.wdl" as RunHISAT2
import "HISAT2RSEMPipeline.wdl" as RunHISAT2RSEM

workflow SmartSeq2SingleCell {
  meta {
    description: "Process SmartSeq2 scRNA-Seq data, include reads alignment, QC metrics collection, and gene expression quantitication"
  }

  # load annotation
  File gtf_file
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat

  # load index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index

  # ref index name
  File hisat2_ref_name
  File hisat2_ref_trans_name

  # samples
  String stranded
  String sample_name
  String output_name
  File fastq1
  File fastq2

  parameter_meta {
    gtf_file: ""
    genome_ref_fasta: ""
    rrna_intervals: ""
    gene_ref_flat: ""
    hisat2_ref_index: ""
    hisat2_ref_trans_index: ""
    rsem_ref_index: ""
    hisat2_ref_name: ""
    hisat2_ref_trans_name: ""
    stranded: ""
    sample_name: ""
    output_name: ""
    fastq1: ""
    fastq2: ""
  }

  call RunHISAT2.RunHISAT2Pipeline as QualityControlWorkflow {
    input:
      fastq_read1 = fastq1,
      fastq_read2 = fastq2, 
      gtf = gtf_file,
      stranded = stranded,
      ref_fasta = genome_ref_fasta,
      rrna_interval = rrna_intervals,
      ref_flat = gene_ref_flat,
      hisat2_ref = hisat2_ref_index,
      hisat2_ref_name = hisat2_ref_name,
      sample_name = sample_name,
      output_prefix = output_name + "_qc"
      
  }

  call RunHISAT2RSEM.RunHISAT2RSEMPipeline as DataWorkflow {
    input:
      fastq_read1 = fastq1, 
      fastq_read2 = fastq2, 
      hisat2_ref_trans = hisat2_ref_trans_index,
      hisat2_ref_trans_name = hisat2_ref_trans_name,
      rsem_genome = rsem_ref_index,
      output_prefix = output_name + "_rsem",
      sample_name = sample_name
  }

  output {
    # quality control outputs
    File aligned_bam = QualityControlWorkflow.aligned_bam
    File alignment_summary_metrics = QualityControlWorkflow.alignment_summary_metrics
    File bait_bias_detail_metrics = QualityControlWorkflow.bait_bias_detail_metrics
    File bait_bias_summary_metrics = QualityControlWorkflow.bait_bias_summary_metrics
    File base_call_dist_metrics = QualityControlWorkflow.base_call_dist_metrics
    File base_call_pdf = QualityControlWorkflow.base_call_pdf
    File dedup_metrics = QualityControlWorkflow.dedup_metrics
    File error_summary_metrics = QualityControlWorkflow.error_summary_metrics
    File gc_bias_detail_metrics = QualityControlWorkflow.gc_bias_detail_metrics
    File gc_bias_dist_pdf = QualityControlWorkflow.gc_bias_dist_pdf
    File gc_bias_summary_metrics = QualityControlWorkflow.gc_bias_summary_metrics
    File insert_size_hist = QualityControlWorkflow.insert_size_hist
    File insert_size_metrics = QualityControlWorkflow.insert_size_metrics
    File hisat2_logfile = QualityControlWorkflow.logfile
    File hisat2_metfile = QualityControlWorkflow.metfile
    File pre_adapter_details_metrics = QualityControlWorkflow.pre_adapter_details_metrics
    File quality_by_cycle_metrics = QualityControlWorkflow.quality_by_cycle_metrics
    File quality_by_cycle_pdf = QualityControlWorkflow.quality_by_cycle_pdf
    File quality_distribution_dist_pdf = QualityControlWorkflow.quality_distribution_dist_pdf
    File quality_distribution_metrics = QualityControlWorkflow.quality_distribution_metrics
    File rna_coverage = QualityControlWorkflow.rna_coverage
    File rna_metrics = QualityControlWorkflow.rna_metrics

    # data outputs
    File aligned_trans_bam = DataWorkflow.aligned_trans_bam
    File hisat2tran_logfile = DataWorkflow.logfile
    File hisat2tran_metfile = DataWorkflow.metfile
    File rsem_cnt_log = DataWorkflow.rsem_cnt_log
    File rsem_gene_results = DataWorkflow.rsem_gene_results
    File rsem_isoform_results = DataWorkflow.rsem_isoform_results
    File rsem_model_log = DataWorkflow.rsem_model_log
    File rsem_theta_log = DataWorkflow.rsem_theta_log
    File rsem_time_log = DataWorkflow.rsem_time_log
  }
}
