task BuildHISAT2forRSEM {
  File rsem_index
  String ref_name
  
  command {
    ##extract rsem index
    tar -xvf ${rsem_index}
    ##building index
    hisat2-build -p 8 rsem/rsem_trans_index.idx.fa  ${ref_name}
    mkdir ${ref_name}
    mv *.ht2 ${ref_name}
    tar -zcvf "${ref_name}.tar.gz" "${ref_name}"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }
  output {
    File hisat2Ref = "${ref_name}.tar.gz"
  }
}
workflow HISAT2Ref {
  String ref_name
  File  rsem_index
  
  call BuildHISAT2forRSEM {
    input:
      ref_name = ref_name,
      rsem_index = rsem_index			
  }
  output {
    File hisat2_ref = BuildHISAT2forRSEM.hisat2Ref
  }
}
