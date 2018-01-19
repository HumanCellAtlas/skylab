import "Kallisto.wdl" as kallisto

workflow kallisto_mkref {
  File transcriptome_fasta  # fasta file containing transcripts ONLY (no genome!)
  Int k  # the size of kmer to build the index using.

  call kallisto.Mkref {
    input:
      transcriptome_fasta = transcriptome_fasta
  }

  output {
    File index = Mkref.index
  }
}
