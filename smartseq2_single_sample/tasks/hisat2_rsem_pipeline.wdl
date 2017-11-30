import "hisat2.wdl" as hisat2
import "rsem.wdl" as rsem

workflow RunHisat2RsemPipeline {
  File fastq_read1
  File fastq_read2
  File hisat2_ref_trans
  File rsem_genome
  String output_prefix
  String hisat2_ref_trans_name
  String sample_name
  
  call hisat2.HISAT2rsem as Hisat2Trans {
    input:
      hisat2_ref = hisat2_ref_trans,
      fq1 = fastq_read1,
      fq2 = fastq_read2,
      ref_name = hisat2_ref_trans_name,
      sample_name = sample_name,
      output_name = output_prefix
    }
  call rsem.RsemExpression as Rsem {
    input:
      trans_aligned_bam = Hisat2Trans.output_bam,
      rsem_genome = rsem_genome,
      rsem_out = output_prefix
    }
  output {
    File aligned_trans_bam = Hisat2Trans.output_bam
    File metfile = Hisat2Trans.metfile
    File logfile = Hisat2Trans.logfile
    File rsem_gene_results = Rsem.rsem_gene
    File rsem_isoform_results = Rsem.rsem_isoform
    File rsem_gene_counts = Rsem.rsem_gene_count
    File rsem_time_log = Rsem.rsem_time
    File rsem_cnt_log = Rsem.rsem_cnt
    File rsem_model_log = Rsem.rsem_model
    File rsem_theta_log = Rsem.rsem_theta
    
  }
}
