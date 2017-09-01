

task align_star {
  File fastq_archive
  File index  # does gcloud require that I specify each file in here?
              # if so, I'll tar it so I can drag the whole thing in.
  File output_directory

  command {
    STAR --runMode alignReads \
      --genomeDir ${index} \
      --outSamUnaligned within \
      --outPrefix ${output_directory} \
  }

  runtime {
    docker: ""
  }

  output {
    File samfile = "${output_directory}/Aligned.out.sam"  # double check this
  }

}


task extract_aligned {

  command {
    ""
  }

  runtime {
    docker: ""
  }

  output {

  }

}


task extract_unaligned {

  command {

  }

  runtime {
    docker: ""
  }

  output {

  }

}


task concatenate {

  command {

  }

  runtime {
    docker: ""
  }

  output {

  }

}


task revert_fastq {

  command {

  }

  runtime {
    docker: ""
  }

  output {

  }

}


workflow generate_test {

  # inputs go here

  call align_star {
    input:
  }

  call extract_unaligned {
    input:
  }

  call extract_aligned {
    input:
  }

  call concatenate {
    input:
  }

  call revert_to_fastq {
    input:
  }

  output {

  }

}