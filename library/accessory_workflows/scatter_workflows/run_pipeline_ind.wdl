import "star_pipeline.wdl" as run_star
import "hisat2_pipeline.wdl" as run_hisat2
import "hisat2_rsem_pipeline.wdl" as run_hisat2_rsem
import "kallisto_pipeline.wdl" as run_kallisto
## this workflow will only run on pipeline based on workflow_name: star,hisat2, hisat2rsem
workflow run_pipelines {

  ## load annotation
  File gtffile
  File genome_ref_fasta
  File rrna_intervals
  File gene_ref_flat
  #load index
  File star_ref_index
  File hisat2_ref_index
  File hisat2_ref_trans_index
  File rsem_ref_index
  File kallisto_ref_index
  ##ref index name
  File hisat2_ref_name
  File hisat2_ref_trans_name
  ## array of string represents featuretypes
  Array[String] featureType = ['exon','gene','transcript']
  ##samples
  String stranded
  File sampleinfo
  String sra_dir
  Array[String] sraIDs = read_lines(sampleinfo)
  ## workflow name
  String workflow_name
  
  scatter(idx in range(length(sraIDs))) { 
    if (workflow_name == "Star" ) {
    call run_star.RunStarPipeline {
      input:
        fastq_read1 = sra_dir+'/'+sraIDs[idx]+"_1.fastq.gz",
        fastq_read2 = sra_dir+'/'+sraIDs[idx]+"_2.fastq.gz",
        gtf = gtffile,
        stranded = stranded,
        ref_fasta = genome_ref_fasta,
        rrna_interval = rrna_intervals,
        ref_flat = gene_ref_flat,
        star_genome = star_ref_index,
        rsem_genome = rsem_ref_index,
        sample_name = sraIDs[idx],
        output_prefix = sraIDs[idx]     
    }
}
  if( workflow_name == "HISAT2" ) { 
   call run_hisat2.RunHisat2Pipeline {
      input:
        fastq_read1 = sra_dir+'/'+sraIDs[idx]+"_1.fastq.gz",
        fastq_read2 = sra_dir+'/'+sraIDs[idx]+"_2.fastq.gz", 
        gtf = gtffile,
        stranded = stranded,
        ref_fasta = genome_ref_fasta,
        rrna_interval = rrna_intervals,
        ref_flat = gene_ref_flat,
        hisat2_ref = hisat2_ref_index,
        hisat2_ref_name = hisat2_ref_name,
        sample_name = sraIDs[idx],
        output_prefix = sraIDs[idx]
    }
}
  if(workflow_name == "HISAT2Rsem") {
    call run_hisat2_rsem.RunHisat2RsemPipeline {
      input:
        fastq_read1 = sra_dir+'/'+sraIDs[idx]+"_1.fastq.gz", 
        fastq_read2 = sra_dir+'/'+sraIDs[idx]+"_2.fastq.gz", 
        hisat2_ref_trans = hisat2_ref_trans_index,
        hisat2_ref_trans_name = hisat2_ref_trans_name,
        rsem_genome = rsem_ref_index,
        output_prefix = sraIDs[idx],
        sample_name = sraIDs[idx]
    }
}
  if(workflow_name == "Kallisto" ) {
    call run_kallisto.RunKallisto {
      input:
        fastq1 = sra_dir+'/'+sraIDs[idx]+"_1.fastq.gz",
        fastq2 = sra_dir+'/'+sraIDs[idx]+"_2.fastq.gz",
        ref_genome = kallisto_ref_index,
        sample_name = sraIDs[idx]
    }
  }
}
}
