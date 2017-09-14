task BuildRsemReference {
  File ref_fasta
  File gtf_file

  command {
    mkdir rsem
    rsem-prepare-reference --gtf ${gtf_file} --bowtie ${ref_fasta} rsem/rsem_trans_index
    tar -cvf rsem.tar rsem/
  }
  output {
    File rsemRef = "rsem.tar"
  }
  runtime {
    docker: "humancellatlas/rsem"
    memory: "10 GB"
    disks: "local-disk 100 HDD"
  }
}

workflow RsemRef {
  File fasta
  File gtf
  
  call BuildRsemReference{
    input:
      ref_fasta = fasta,
      gtf_file = gtf
  }
  output {
    File rsem_ref = BuildRsemReference.rsemRef
  }
}
