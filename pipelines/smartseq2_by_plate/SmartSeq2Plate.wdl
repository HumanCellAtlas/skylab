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
task AggregatePicardMetricsByCell{
  Array[File] metric_files
  String output_name
  String docker = "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
  command {
    set -e
    gsutil cp gs://broad-dsde-mint-dev-teststorage/pipeline_testing_scripts/AggregatePicardMetricsByCell.py ./
    python AggregatePicardMetricsByCell.py -f ${sep=',' metric_files}   -o ${output_name}
  }
  output{
    File aggregated_result = output_name+".picard.core.json"
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

task AggregatePicardMetricsByBatch {
  Array[File] metric_files
  String output_name
  String docker = "gcr.io/broad-dsde-mint-dev/benchmarking-tools:0.0.1"
  command {
    set -e
    gsutil cp gs://broad-dsde-mint-dev-teststorage/pipeline_testing_scripts/AggregatePicardMetrics.py ./
    python AggregatePicardMetrics.py -f ${sep=',' metric_files}   -o ${output_name}
  }
  output{
    File aggregated_result = output_name+".picard.batch.core.csv"
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
   call AggregatePicardMetricsByCell {
      input:
        metric_files = [sc.alignment_summary_metrics,sc.rna_metrics, sc.dedup_metrics,sc.insert_size_metrics],
        output_name = sraIDs[idx]
      } 
  }
  scatter(i in range(length(out_types))){
    call AggregateDataMatrix as AggregateGene {
      input:
        filename_array = sc.rsem_gene_results,
        col_name = out_types[i],
        output_name ="gene_" + out_types[i] + ".csv"
      }
    call AggregateDataMatrix as AggregateIsoform {
      input:
        filename_array = sc.rsem_isoform_results,
        col_name = out_types[i],
        output_name = "isoform_" + out_types[i] + ".csv"
      }
    }
  call AggregatePicardMetricsByBatch {
    input:
      metric_files = AggregatePicardMetricsByCell.aggregated_result,
      output_name = batch_id  
}
}  
