import "SmartSeq2SingleSample.wdl" as single_cell_run

task AggregateDataMatrix{
  Array[File] filename_array
  String col_name
  String output_name
  String docker = "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
  command {
    set -e
    gsutil cp gs://broad-dsde-mint-dev-teststorage/pipeline_testing_scripts/MergeDataMatrix.py ./
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
  String docker = "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
  String run_type
  command {
    set -e
    gsutil cp gs://broad-dsde-mint-dev-teststorage/pipeline_testing_scripts/AggregateMetrics.py ./
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
        output_name =batch_id+"_gene_" + out_types[i] + ".csv"
      }
    call AggregateDataMatrix as AggregateIsoform {
      input:
        filename_array = sc.rsem_isoform_results,
        col_name = out_types[i],
        output_name = batch_id+"_isoform_" + out_types[i] + ".csv"
    }
 }
  Array[String] row_metrics_names = ["alignment_summary_metrics","rna_metrics","dedup_metrics","insert_size_metrics","gc_bias_summary_metrics"]
  Array[String] table_metrics_names = ["base_call_dist_metrics","gc_bias_detail_metrics","pre_adapter_details_metrics","bait_bias_detail_metrics","bait_bias_summary_metrics","error_summary_metrics"]
  Array[Array[File]] row_metrics = [sc.alignment_summary_metrics,sc.rna_metrics, sc.dedup_metrics, sc.insert_size_metrics,sc.gc_bias_summary_metrics]
  Array[Array[File]] table_metrics = [sc.base_call_dist_metrics,sc.gc_bias_detail_metrics,sc.pre_adapter_details_metrics,sc.bait_bias_detail_metrics,sc.bait_bias_summary_metrics,sc.error_summary_metrics]
  scatter(i in range(length(row_metrics))){
    
    call AggregateQCMetrics as AggregateRow {
      input:
        metric_files = row_metrics[i],
        output_name = batch_id+"_"+row_metrics_names[i],
        run_type = "Picard"
    }
  }
  scatter(i in range(length(table_metrics))){
    call AggregateQCMetrics as AppendTable {
      input:
        metric_files = table_metrics[i],
        output_name = batch_id+"_"+table_metrics_names[i],
        run_type = "PicardTable"
    }
  }
 
  call AggregateQCMetrics as AggregateAll {
    input:
      metric_files = AggregateRow.aggregated_result,
      output_name = batch_id+"_"+"QCs",
      run_type = "ALL"
  }
}  
