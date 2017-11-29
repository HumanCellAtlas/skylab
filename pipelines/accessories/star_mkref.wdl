
import "StarMkref.wdl" as star

workflow star_mkref {
  File fasta_file  # genome fasta file
  File annotation_file  # matching .gtf annotation file
  
  call star.StarMkref {
    input:
      fasta_file = fasta_file,
      annotation_file = annotation_file
  }

  output {
    File genome = StarMkref.genome
  }
}
