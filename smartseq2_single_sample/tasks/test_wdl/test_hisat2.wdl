import "hisat2.wdl" as hisat2

workflow test_hisat2 {
  File hisat2_ref
  File fq1
  File fq2
  String ref_name
  String output_name

  call hisat2.HISAT2 as test_hisat2 {
    input:
      hisat2_ref = hisat2_ref,
      fq1 = fq1,
      fq2 = fq2,
      ref_name = ref_name,
      output_name = output_name 
    }
  output {
    File logFile = test_hisat2.logfile
    File metfile = test_hisat2.metfile
    File bamfile = test_hisat2.output_bam
  }
}
