import "hisat2.wdl" as hisat2

workflow test_hisat2 {
  File hisat2_ref
  File fq1
  File fq2
  String ref_name
  String output_name
  String sample_name
 ## to inspect reference build, such as reference, num of exon, num of snps and num of splicing
  call hisat2.hisat2_inspect_index as test_hisat2_inspect {
    input:
      hisat2_ref = hisat2_ref,
      ref_name = ref_name
  }
  ## run paired-end reads with HISAT2
  call hisat2.HISAT2PE as test_hisat2_pe {
    input:
      hisat2_ref = hisat2_ref,
      fq1 = fq1,
      fq2 = fq2,
      ref_name = ref_name,
      sample_name =sample_name,
      output_name = output_name
  }
  ## Run single end reads with HISAT2
  call hisat2.HISAT2SE as test_hisat2_se {
    input:
      hisat2_ref = hisat2_ref,
      fq = fq1,
      ref_name = ref_name,
      sample_name = sample_name,
      output_name = output_name,
   }
  output {
    File logFile_pe = test_hisat2_pe.logfile
    File metfile_pe = test_hisat2_pe.metfile
    File bamfile_pe = test_hisat2_pe.output_bam
    File logFile_se = test_hisat2_se.logfile
    File metfile_se = test_hisat2_se.metfile
    File bamfile_se = test_hisat2_se.output_bam
    File inspectlog = test_hisat2_inspect.logfile
  }
}
