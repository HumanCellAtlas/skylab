task BuildHISAT2reference{
  String ref_name
  String gtf_version
  String dbsnp_version

  command {
      
    /opt/tools/hisat2-2.1.0/make_grch38_snp_tran_ensembl.sh 90 150 ${ref_name}
  }
  output {
    File hisat2Ref = "${ref_name}.tar.gz"
  }
  runtime {
    docker:"humancellatlas/hisat2:2-2.1.0"
    memory: "250 GB"
    disks: "local-disk 50 HDD"
    cpu: "8"
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
