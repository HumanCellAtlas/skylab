version 1.0

task BuildStarReference {
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
  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"

  command <<<
    set -eo pipefail

    ## download fasta
    wget ~{ftp_path}/~{genome_fa}.gz
    gunzip ~{genome_fa}.gz

    ## download gtf file
    wget ~{ftp_path}/~{annotation_gtf}.gz
    gunzip ~{annotation_gtf}.gz

    mkdir star
    STAR --runMode genomeGenerate \
      --genomeDir star \
      --genomeFastaFiles ~{genome_fa} \
      --sjdbGTFfile ~{annotation_gtf} \
      --sjdbOverhang 100 \
      --runThreadN 16
    tar -cvf "~{ref_name}.tar" star
  >>>

  output {
    File starRef = "~{ref_name}.tar"
  }

  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-1.0.0"
    memory: "50 GB"
    disks :"local-disk 100 HDD"
    cpu:"16"
  }
}

workflow StarRef {
  input {
    String gtf_version
    String organism
    String organism_prefix
  }

  call BuildStarReference {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix
  }

  output {
    File star_ref = BuildStarReference.starRef
  }
}
