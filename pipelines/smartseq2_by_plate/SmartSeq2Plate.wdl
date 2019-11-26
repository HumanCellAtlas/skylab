import "SmartSeq2SingleSample.wdl" as single_cell_run
import "SmartSeq2PlateAggregation.wdl" as ss2_plate_aggregation
       
workflow RunSmartSeq2ByPlate {
  # load annotation
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat
  File gtf_file

  # load index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index

  # ref index name
  File hisat2_ref_name
  File hisat2_ref_trans_name

  # samples
  String stranded
  String sra_dir
  Array[String] out_types = ["est_counts","tpm"]  
  Array[String] sraIDs
  String batch_id
  String docker

  # Run the SS2 pipeline for each cell in a scatter
  scatter(idx in range(length(sraIDs))) { 
    call single_cell_run.SmartSeq2SingleCell as sc {
      input:
        fastq1 = sra_dir + '/' + sraIDs[idx] + "_1.fastq.gz",
        fastq2 = sra_dir + '/' + sraIDs[idx] + "_2.fastq.gz", 
        gtf_file = gtf_file,
        stranded = stranded,
        genome_ref_fasta = genome_ref_fasta,
        rrna_intervals = rrna_intervals,
        gene_ref_flat = gene_ref_flat,
        hisat2_ref_index = hisat2_ref_index,
        hisat2_ref_name = hisat2_ref_name,
        hisat2_ref_trans_index = hisat2_ref_trans_index,
        hisat2_ref_trans_name = hisat2_ref_trans_name,
        rsem_ref_index = rsem_ref_index,
        sample_name = sraIDs[idx],
        output_name = sraIDs[idx]
   }
  }

  call ss2_plate_aggregation.AggregateDataMatrix as AggregateGene {
    input:
      filename_array = sc.rsem_gene_results,
        col_name = "est_counts",
        output_name =batch_id+"_gene_est_counts.csv",
        docker = docker
  }

   call ss2_plate_aggregation.AggregateDataMatrix as AggregateIsoform {
      input:
        filename_array = sc.rsem_isoform_results,
        col_name = "tpm",
        output_name = batch_id+"_isoform_tpm.csv",
        docker = docker
    }
 
  Array[String] row_metrics_names = ["alignment_summary_metrics","rna_metrics","dedup_metrics","insert_size_metrics","gc_bias_summary_metrics"]
  
  Array[String] table_metrics_names = ["base_call_dist_metrics","gc_bias_detail_metrics","pre_adapter_details_metrics","bait_bias_detail_metrics","bait_bias_summary_metrics","error_summary_metrics"]
  
  Array[Array[File]] row_metrics = [sc.alignment_summary_metrics,sc.rna_metrics, sc.dedup_metrics, sc.insert_size_metrics,sc.gc_bias_summary_metrics]
  
  Array[Array[File]] table_metrics = [sc.base_call_dist_metrics,sc.gc_bias_detail_metrics,sc.pre_adapter_details_metrics,sc.bait_bias_detail_metrics,sc.bait_bias_summary_metrics,sc.error_summary_metrics]

  Array[String] hisat2_logs_names = ["hisat2_log_file","hisat2_transcriptome_log_file"]
  Array[Array[File]] hisat2_logs = [sc.hisat2_log_file,sc.hisat2_transcriptome_log_file]

  Array[String] rsem_logs_names = ["rsem_cnt_log"] 
  Array[Array[File]] rsem_logs = [sc.rsem_cnt_log]

  scatter(i in range(length(row_metrics))){
    call ss2_plate_aggregation.AggregateQCMetrics as AggregateRow {
      input:
        metric_files = row_metrics[i],
        output_name = batch_id+"_"+row_metrics_names[i],
        run_type = "Picard",
        docker = docker
    }
  }
  
  scatter(i in range(length(hisat2_logs))){
    call ss2_plate_aggregation.AggregateQCMetrics as AggregateHisat2 {
      input:
        metric_files = hisat2_logs[i],
        output_name = batch_id+"_"+hisat2_logs_names[i],
        run_type = "HISAT2",
        docker  = docker
    }
  }

 scatter(i in range(length(rsem_logs))){
    call ss2_plate_aggregation.AggregateQCMetrics as AggregateRsem {
      input:
        metric_files = rsem_logs[i],
        output_name = batch_id+"_"+rsem_logs_names[i],
        run_type = "RSEM",
        docker = docker
    }
  } 

  scatter(i in range(length(table_metrics))){
    call ss2_plate_aggregation.AggregateQCMetrics as AppendTable {
      input:
        metric_files = table_metrics[i],
        output_name = batch_id+"_"+table_metrics_names[i],
        run_type = "PicardTable",
        docker = docker
    }
  }
 
  call ss2_plate_aggregation.AggregateQCMetricsCore as AggregateCore {
    input:
      picard_metric_files = AggregateRow.aggregated_result,
      hisat2_stats_files= AggregateHisat2.aggregated_result,
      rsem_stats_files = AggregateRsem.aggregated_result,
      output_name = batch_id+"_"+"QCs",
      run_type = "Core",
      docker = docker
  }

  output {
    File core_QC = AggregateCore.aggregated_result
    Array[File] qc_tabls = AppendTable.aggregated_result
    Array[File] gene_matrix = AggregateGene.aggregated_result
    Array[File] isoform_matrix = AggregateIsoform.aggregated_result
  }
}  
