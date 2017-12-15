import "Kallisto.wdl" as kallisto

workflow RunKallisto {
  File fastq1
  File fastq2
  File ref_genome
  String sample_name

  call kallisto.QuantPairedEndNoBam {
    input:
      r1 = fastq1,
      r2 = fastq2,
      index = ref_genome,
      sample_name = sample_name
  }
  output {
    File count_h5 = QuantPairedEndNoBam.abundance_h5
    File count_tsv = QuantPairedEndNoBam.abundance_tsv
    File logfile = QuantPairedEndNoBam.log
  }
}
