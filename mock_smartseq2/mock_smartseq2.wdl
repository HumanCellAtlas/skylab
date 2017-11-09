task calc_expression {
  File left_fastq
  File right_fastq
  File ref_tar

  command {
    echo "Aligning fastqs and calculating expression"
    tar xf ${ref_tar}
    rsem-calculate-expression -p 1 --paired-end ${left_fastq} ${right_fastq} hg19_ucsc_genomeStudio_genes/rsem_trans_index output
    cut -f 1,5 output.isoforms.results > output_matrix.txt
    echo "Created output.transcript.bam and output_matrix.txt"
    echo "Creating qc matrix"
    cp output_matrix.txt qc_matrix.txt
    echo "Wrote qc matrix"
  }
  runtime {
    docker: "humancellatlas/rsem:v1.3.0"
    memory: "2 GB"
    disks: "local-disk 10 HDD"
  }
  output {
    File bam = "output.transcript.bam"
    File matrix = "output_matrix.txt"
    File qc_matrix = "qc_matrix.txt"
  }
}

workflow mock_smartseq2 {
  File left_fastq
  File right_fastq
  File ref_tar

  call calc_expression {
    input:
      left_fastq = left_fastq,
      right_fastq = right_fastq,
      ref_tar = ref_tar
  }
}
