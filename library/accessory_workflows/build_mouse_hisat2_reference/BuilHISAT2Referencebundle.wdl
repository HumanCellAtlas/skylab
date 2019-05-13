task BuildHISAT2reference{
  String ref_name ## name of the tar.gz file without tar.gz suffix
  String gtf_version  ## the actually number of gencode, ex.  27
  String dbsnp_version ## dbsnp version, integer num, ex 150
  
  command {
    make_grch38_snp_tran_gencode.sh ${gtf_version} ${dbsnp_version} 
    mkdir ${ref_name}
    cp *.ht2 ${ref_name}
    tar -zcvf "${ref_name}.tar.gz" "${ref_name}"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "240 GB"
    disks: "local-disk 100 HDD"
    cpu: "16"
  }
  output {
    File hisat2Ref = "${ref_name}.tar.gz"
  }
}

workflow HISAT2Ref {
  String ref_name
  String gtf_version
  String dbsnp_version

  call BuildHISAT2reference {
    input:
      ref_name = ref_name,
      gtf_version = gtf_version,
      dbsnp_version = dbsnp_version
  }
  output {
    File hisat2_ref = BuildHISAT2reference.hisat2Ref
  }
}
