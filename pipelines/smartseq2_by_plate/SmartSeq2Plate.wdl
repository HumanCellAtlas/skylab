import "SmartSeq2SingleSample.wdl" as single_cell_run

task AggregateDataMatrix{
  Array[File] filename_array
  String col_name
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.0.2"
  command {
    set -e
    git clone --branch jx-ss2-platebundle https://github.com/HumanCellAtlas/skylab
    cd skylab/pipelines/smartseq2_by_plate/
    python MergeDataMatrix.py -f ${sep=',' filename_array}  -t ${col_name} -o ${output_name}
  }
  output{
    File aggregated_result = "${output_name}"
  }
  runtime {
    docker: docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 5
    maxRetries: 1
  }
}
task AggregateQCMetrics{
  Array[File] metric_files
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.0.2"
  String run_type
  command {
    set -e
    git clone --branch jx-ss2-platebundle https://github.com/HumanCellAtlas/skylab
    cd skylab/pipelines/smartseq2_by_plate/
    python AggregateMetrics.py -f ${sep=' ' metric_files}   -o ${output_name} -t ${run_type}
  }
  output{
    File aggregated_result = output_name+".csv"
 }
  runtime {
    docker: docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 5
    maxRetries: 1
  }
}
task AggregateQCMetricsCore{
  Array[File] picard_metric_files
  Array[File] rsem_stats_files
  Array[File] hisat2_stats_files
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-python3-scientific:0.0.2"
  String run_type
  command {
    set -e
    git clone --branch jx-ss2-platebundle https://github.com/HumanCellAtlas/skylab
    cd skylab/pipelines/smartseq2_by_plate/
    python AggregateMetrics.py -f ${sep=' ' picard_metric_files} ${sep=' ' hisat2_stats_files} ${sep=' ' rsem_stats_files}   -o ${output_name} -t ${run_type}
  }
  output{
    File aggregated_result = output_name+".csv"
 }
  runtime {
    docker: docker
    memory: "2 GB"
    disks: "local-disk 10 HDD"
    cpu: 1
    preemptible: 5
    maxRetries: 1
  }
}


task MergeBamFiles {
  Array[File] bam_files
  String output_name
  String docker = "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
  
  command {
    set -e 
    samtools merge -@ 8 "${output_name}.bam" ${sep=' ' bam_files} 
    samtools index "${output_name}.bam"
    ls .    
  }
  output {
    File output_bam = output_name + ".bam"
    File output_index = output_name + ".bam.bai"
  }
  runtime {
    docker: docker
    memory: "15 GB"
    disks: "local-disk 100 HDD"
    cpu: 8
    preemptible: 5
    maxRetries: 1
  }
}

workflow RunSmartSeq2ByPlate {

  # load annotation
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat
  File gtf_file
  #load index
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
  scatter(i in range(length(out_types))){
    call AggregateDataMatrix as AggregateGene {
      input:
        filename_array = sc.rsem_gene_results,
        col_name = out_types[i],
        output_name =batch_id+"_gene_" + out_types[i] + ".csv",
        docker = docker
      }
    call AggregateDataMatrix as AggregateIsoform {
      input:
        filename_array = sc.rsem_isoform_results,
        col_name = out_types[i],
        output_name = batch_id+"_isoform_" + out_types[i] + ".csv",
        docker = docker
    }
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
    
    call AggregateQCMetrics as AggregateRow {
      input:
        metric_files = row_metrics[i],
        output_name = batch_id+"_"+row_metrics_names[i],
        run_type = "Picard",
        docker = docker
    }
  }
  scatter(i in range(length(hisat2_logs))){
    
    call AggregateQCMetrics as AggregateHisat2 {
      input:
        metric_files = hisat2_logs[i],
        output_name = batch_id+"_"+hisat2_logs_names[i],
        run_type = "HISAT2",
        docker  = docker
    }
  }

 scatter(i in range(length(rsem_logs))){
  
    call AggregateQCMetrics as AggregateRsem {
      input:
        metric_files = rsem_logs[i],
        output_name = batch_id+"_"+rsem_logs_names[i],
        run_type = "RSEM",
        docker = docker
    }
  } 
  scatter(i in range(length(table_metrics))){
    call AggregateQCMetrics as AppendTable {
      input:
        metric_files = table_metrics[i],
        output_name = batch_id+"_"+table_metrics_names[i],
        run_type = "PicardTable",
        docker = docker
    }
  }
 
  call AggregateQCMetricsCore as AggregateCore {
    input:
      picard_metric_files = AggregateRow.aggregated_result,
      hisat2_stats_files= AggregateHisat2.aggregated_result,
      rsem_stats_files = AggregateRsem.aggregated_result,
      output_name = batch_id+"_"+"QCs",
      run_type = "Core",
      docker = docker
  }
 call MergeBamFiles {
    input:
      bam_files = sc.aligned_bam,
      output_name = batch_id
 }
 output {
  File merged_bam = MergeBamFiles.output_bam
  File bam_index = MergeBamFiles.output_index
  File core_QC = AggregateCore.aggregated_result
  Array[File] qc_tabls = AppendTable.aggregated_result
  Array[File] gene_matrix = AggregateGene.aggregated_result
  Array[File] isoform_matrix = AggregateIsoform.aggregated_result
  }
}  
