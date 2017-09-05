import "ss2_single_sample_wf.wdl" as singlesample
task GatherMetricsToTable {
  Array[File] input_metrics_list
  String output_filename

  command {
    python /tools/scripts/parse_picard.py -I ${write_lines(input_metrics_list)} -T "table" -O "${output_filename}"
  }
  runtime {
    docker: "humancellatlas/python3-application:v0.1"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
  output {
    File metrics_table = "${output_filename}_metrics.csv"
  }
}

workflow Ss2RunMultiSample {
  File sra_list_file
  File gtf
  File ref_fasta
  File rrna_interval
  File ref_flat
  String star_genome
  String rsem_genome
  String sra_dir
  
## start to scatter single sample workflow by sraID
  Array[String] sraIDs=read_lines(sra_list_file)
   
  scatter(idx in range(length(sraIDs))) {
    call singlesample.Ss2RunSingleSample as single_run {  
      input:
        fastq_read1 = sra_dir+'/'+sraIDs[idx]+"_1.fastq.gz",
        fastq_read2 = sra_dir+'/'+sraIDs[idx]+"_2.fastq.gz",
        gtf = gtf,
        ref_fasta = ref_fasta,
        rrna_interval = rrna_interval,
        ref_flat = ref_flat,
        star_genome = star_genome,
        rsem_genome = rsem_genome,
        output_prefix = sraIDs[idx]
    
    }
 }

  call GatherMetricsToTable as collect_all {
    input:
      input_metrics_list = single_run.sample_metrics,
      output_filename = "all_samples"
  }
  output {
    single_run.*
    collect_all.*
  }
}
