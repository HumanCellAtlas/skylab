version 1.0

task BuildHISAT2Trans {
  input {
    String gtf_version
    # Either 'human' or 'mouse'
    String organism
    # Either 'h' or 'm'
    String organism_prefix
  }

  String ftp_path = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_~{organism}/release_~{gtf_version}"
  String genome_fa = "GRC~{organism_prefix}38.primary_assembly.genome.fa"
  String annotation_gtf = "gencode.v~{gtf_version}.primary_assembly.annotation.gtf"
  String ref_name = "gencode_~{organism}_v~{gtf_version}_trans"

  command <<<
    set -eo pipefail

    ## download fasta
    wget ~{ftp_path}/~{genome_fa}.gz
    gunzip ~{genome_fa}.gz

    ## download gtf file
    wget ~{ftp_path}/~{annotation_gtf}.gz
    gunzip ~{annotation_gtf}.gz

    ## extract transcriptome fa
    gffread -w transcriptome.fa -g ~{genome_fa} ~{annotation_gtf}

    ## building index
    hisat2-build -p 8 transcriptome.fa ~{ref_name}
    mkdir ~{ref_name}
    mv *.ht2 ~{ref_name}
    tar -zcvf ~{ref_name}.tar.gz ~{ref_name}
  >>>

  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }
  output {
    File hisat2Ref = "~{ref_name}.tar.gz"
  }
}

workflow HISAT2Ref {
  input {
    String gtf_version
    String organism
    String organism_prefix
   }
  
  call BuildHISAT2Trans {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix
  }
  output {
    File hisat2_ref = BuildHISAT2Trans.hisat2Ref
  }
}
