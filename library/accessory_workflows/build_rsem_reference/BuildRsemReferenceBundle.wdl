task BuildRsemReference {
  File ref_fasta
  File gtf_file
  String ref_name
  command {
    mkdir rsem
    rsem-prepare-reference --gtf ${gtf_file} --bowtie ${ref_fasta} rsem/rsem_trans_index
    tar -cvf ${ref_name}.tar rsem/
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0"
    memory: "10 GB"
    disks: "local-disk 100 HDD"
  }
  output {
    File rsemRef = "${ref_name}.tar"
  }
}

workflow RsemRef {
  File fasta
  File gtf
  String ref_name

  call BuildRsemReference {
    input:
      ref_fasta = fasta,
      gtf_file = gtf,
      ref_name = ref_name
  }
  output {
    File rsem_ref = BuildRsemReference.rsemRef
  }
}
