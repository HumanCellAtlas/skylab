##### Naming Rules#####
## ** Files should be named in CamelCase starting with a capital
## ** Tasks should be named in CamelCase starting with a capital
## ** Alias' should be named in CamelCase starting with a capital
## ** Variable declarations should be lower_snake_case

##### New Line Rules #####
## ** New line before and after `calls`
## ** New lines between blocks in tasks
## ** New line variable declarations that you want to group
## ** New line between imports and workflow definition
## ** No new lines at the start or end of any workflow/call/task definition
## ** No need for having two new lines anywhere
## ** No new lines at the start of any blocks

import "Hisat2QcPipeline.wdl" as RunHisatQCWDL  ## ** Should import alias follow some kind of pattern like ending in "WDL" or something else
import "Hisat2RsemPipeline.wdl" as RunHisatRSEMWDL

## ** Explanations above the workflow about what inputs it expects, what outputs/metrics it outputs
## ** What, if any, assumptions are made in the pipeline
## ** In general any kind of description you want to add above the workflow is nice
## ** A good example of this can be found in these wdls
## ** https://github.com/broadinstitute/dsde-pipelines/blob/develop/genomes_in_the_cloud/single_sample/PairedSingleSampleWf.wdl
## ** https://github.com/broadinstitute/gatk/blob/master/scripts/mutect2_wdl/mutect2.wdl

## main pipeline:
## QC track: HISAT2+Picard
## this track will produce aligned bam file and a set of QC metrics
## rsem quantification track: HISAT2+RSEM
## this track involves hisat2 transcriptome alignment and RSEM gene expression estimation.
workflow SmartSeq2SingleCell {
  # load annotation
  File gtf_file
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat
  #load index
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

  call RunHisatQCWDL.RunHisat2Pipeline as QC {
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

  call RunHisatRSEMWDL.RunHisat2RsemPipeline as Data {
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
    ## qc outputs
    File aligned_bam = QC.aligned_bam
    File alignment_summary_metrics = QC.alignment_summary_metrics
    File bait_bias_detail_metrics = QC.bait_bias_detail_metrics
    File bait_bias_summary_metrics = QC.bait_bias_summary_metrics
    File base_call_dist_metrics = QC.base_call_dist_metrics
    File base_call_pdf = QC.base_call_pdf
    File dedup_metrics = QC.dedup_metrics
    File error_summary_metrics = QC.error_summary_metrics
    File gc_bias_detail_metrics = QC.gc_bias_detail_metrics
    File gc_bias_dist_pdf = QC.gc_bias_dist_pdf
    File gc_bias_summary_metrics = QC.gc_bias_summary_metrics
    File insert_size_hist = QC.insert_size_hist
    File insert_size_metrics = QC.insert_size_metrics
    File hisat2_logfile = QC.logfile
    File hisat2_metfile = QC.metfile
    File pre_adapter_details_metrics = QC.pre_adapter_details_metrics
    File quality_by_cycle_metrics = QC.quality_by_cycle_metrics
    File quality_by_cycle_pdf = QC.quality_by_cycle_pdf
    File quality_distribution_dist_pdf = QC.quality_distribution_dist_pdf
    File quality_distribution_metrics = QC.quality_distribution_metrics
    File rna_coverage = QC.rna_coverage
    File rna_metrics = QC.rna_metrics
    ## data outputs
    File aligned_trans_bam = Data.aligned_trans_bam
    File hisat2tran_logfile = Data.logfile
    File hisat2tran_metfile = Data.metfile
    File rsem_cnt_log = Data.rsem_cnt_log
    File rsem_gene_results = Data.rsem_gene_results
    File rsem_isoform_results = Data.rsem_isoform_results
    File rsem_model_log = Data.rsem_model_log
    File rsem_theta_log = Data.rsem_theta_log
    File rsem_time_log = Data.rsem_time_log
  }
}
