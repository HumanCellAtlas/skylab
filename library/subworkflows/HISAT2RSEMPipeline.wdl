import "HISAT2.wdl" as HISAT2
import "RSEM.wdl" as RSEM

workflow RunHISAT2RSEMPipeline {
  meta {
    description: "JISHU HELP AGAINNNNNNNNNNNNN"
  }

  File fastq_read1
  File fastq_read2
  File hisat2_ref_trans
  File rsem_genome
  String output_prefix
  String hisat2_ref_trans_name
  String sample_name

  parameter_meta {
    fastq_read1: ""
    fastq_read2: ""
    hisat2_ref_trans: ""
    rsem_genome: ""
    output_prefix: ""
    hisat2_ref_trans_name: ""
    sample_name: ""
  }
  ################### TODO should `Trans` be Transcript?? ##################
  call HISAT2.HISAT2RSEM as HISAT2Trans {
    input:
      hisat2_ref = hisat2_ref_trans,
      fq1 = fastq_read1,
      fq2 = fastq_read2,
      ref_name = hisat2_ref_trans_name,
      sample_name = sample_name,
      output_name = output_prefix,
  }

  call RSEM.RSEMExpression {
    input:
      trans_aligned_bam = Hisat2Trans.output_bam,
      rsem_genome = rsem_genome,
      rsem_out = output_prefix,
      disk_size = ceil(fastq_size * bam_disk_multiplier+rsem_ref_size+additional_disk * 2.0)
  }

  output {
    File aligned_trans_bam = HISAT2Trans.output_bam
    File metfile = HISAT2Trans.metfile
    File logfile = HISAT2Trans.logfile
    File rsem_gene_results = RSEMExpression.rsem_gene
    File rsem_isoform_results = RSEMExpression.rsem_isoform
    File rsem_time_log = RSEMExpression.rsem_time
    File rsem_cnt_log = RSEMExpression.rsem_cnt
    File rsem_model_log = RSEMExpression.rsem_model
    File rsem_theta_log = RSEMExpression.rsem_theta
  }
}
