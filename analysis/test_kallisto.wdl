import "Kallisto.wdl" as kallisto

workflow test_kallisto {
  File transcriptome_fasta
  File r1
  File r2
  File umi
  Int k

  call kallisto.Mkref {
    input:
      k = k,
      transcriptome_fasta = transcriptome_fasta
  }

  call kallisto.PseudoSingleEndUMI {
    input:
      index = Mkref.index,
      r1 = r1,
      umi = umi
  }

  call kallisto.QuantSingleEnd {
    input:
      index = Mkref.index,
      r1 = r1
  }

  call kallisto.QuantPairedEnd {
    input:
      index = Mkref.index,
      r1 = r1,
      r2 = r2
  }
}