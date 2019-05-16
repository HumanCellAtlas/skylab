version 1.0

task BuildStar {
  input {
    String gtf_version
    String organism
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
    File star_index = "~{ref_name}.tar"
  }

  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-1.0.0"
    memory: "50 GB"
    disks :"local-disk 100 HDD"
    cpu:"16"
  }
}

task BuildRsem {
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

  command {
    set -eo pipefail

    ## download fasta
    wget ~{ftp_path}/~{genome_fa}.gz
    gunzip ~{genome_fa}.gz

    ## download gtf file
    wget ~{ftp_path}/~{annotation_gtf}.gz
    gunzip ~{annotation_gtf}.gz

    mkdir rsem
    rsem-prepare-reference --gtf ~{annotation_gtf} --bowtie ~{genome_fa} rsem/rsem_trans_index
    tar -cvf ~{ref_name}.tar rsem/
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0"
    memory: "10 GB"
    disks: "local-disk 100 HDD"
  }
  output {
    File rsem_index = "${ref_name}.tar"
  }
}

task BuildHisat2FromRsem {
  input {
    File rsem_index
  }

  String rsem_reference_name = basename(rsem_index, ".tar")
  String ref_name = "hisat2_from_rsem_~{rsem_reference_name}"

  command {

    # extract rsem index
    tar -xf ~{rsem_index}

    # build index
    hisat2-build -p 8 rsem/rsem_trans_index.idx.fa "~{ref_name}"
    mkdir "~{ref_name}"
    mv ./*.ht2 "~{ref_name}"
    tar -zcvf "~{ref_name}.tar.gz" "~{ref_name}"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }
  output {
    File hisat2_index = "~{ref_name}.tar.gz"
  }
}

task BuildHisat2 {
  # This method builds a HISAT2 index directly from the downloaded GTF and Fasta files on GENCODE
  # It does not include dbsnp, whereas the below method does.

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

  command {
    set -eo pipefail

    # download fasta
    wget ~{ftp_path}/~{genome_fa}.gz

    # download gtf file
    wget ~{ftp_path}/~{annotation_gtf}.gz
    gunzip ~{annotation_gtf}.gz

    # build index
    hisat2-build -p 8 ~{genome_fa} ~{ref_name}
    mkdir ~{ref_name}
    mv ./*.ht2 ~{ref_name}
    tar -zcvf "~{ref_name}.tar.gz" "~{ref_name}"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }
  output {
    File hisat2_index = "${ref_name}.tar.gz"
  }
}

task BuildHisat2SnpHaplotypeSplicing {
  # This version includes SNP, haplotype, and
  input {
    String organism  # Either 'human' or 'mouse'
    String organism_prefix  # Either 'h' or 'm'
    String genome_short_string  # e.g. hg38, mm10

    String gtf_version  # the actually number of gencode, ex.  27
    String dbsnp_version  # dbsnp version, integer num, ex 150
  }

  String ftp_path = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_~{organism}/release_~{gtf_version}"
  String genome_fa = "GRC~{organism_prefix}38.primary_assembly.genome.fa"
  String annotation_gtf = "gencode.v~{gtf_version}.primary_assembly.annotation.gtf"
  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String snp_file = "snp~{dbsnp_version}Common.txt"

  command <<<

    HISAT2_DIR=/opt/tools/hisat2-2.1.0

    # download fasta
    wget ~{ftp_path}/~{genome_fa}.gz

    # download gtf file
    wget ~{ftp_path}/~{annotation_gtf}.gz
    gunzip ~{annotation_gtf}.gz

    # download snp file
    wget http://hgdownload.cse.ucsc.edu/goldenPath/~{genome_short_string}/database/"$SNP_FILE".gz
    gunzip ~{snp_file}.gz

    # extract snps, splice sites, and exon information
    $HISAT2_DIR/hisat2_extract_snps_UCSC.py ~{genome_fa} "$SNP_FILE genome"
    $HISAT2_DIR/hisat2_extract_splice_sites.py ~{annotation_gtf} > genome.ss
    $HISAT2_DIR/hisat2_extract_exons.py ~{annotation_gtf} > genome.exon

    # build the hisat2 reference
    $HISAT2_DIR/hisat2-build \
      -p 8 \
      genome.fa \
      --snp genome.snp \
      --haplotype genome.haplotype \
      --ss genome.ss \
      --exon genome.exon \
      genome_snp_tran

    mkdir ~{ref_name}
    cp ./*.ht2 ~{ref_name}
    tar -zcvf "~{ref_name}.tar.gz" "~{ref_name}"

  >>>

  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.3.0-2-2.1.0"
    memory: "240 GB"
    disks: "local-disk 100 HDD"
    cpu: "16"
  }
  output {
    File hisat2_index = "~{ref_name}.tar.gz"
  }
}

task BuildPicardRefFlat {
  input {
    String gtf_version
  }

  command {
    set -eo pipefail

    wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_~{gtf_version}/gencode.v~{gtf_version}.primary_assembly.annotation.gtf.gz
  }

}

workflow BuildIndices {
  input {
    String gtf_version
    String organism
    String organism_prefix
    String genome_short_string
    String dbsnp_version
  }

  call BuildStar {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix
  }

  call BuildRsem {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix
  }

  call BuildHisat2FromRsem {
    input:
      rsem_index = BuildRsem.rsem_index
  }

  call BuildHisat2 {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix
  }

  call BuildHisat2SnpHaplotypeSplicing {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix,
      genome_short_string = genome_short_string,
      dbsnp_version = dbsnp_version
  }
  output {
    File star_index = BuildStar.star_index
    File rsem_index = BuildRsem.rsem_index
    File hisat2_from_rsem_index = BuildHisat2FromRsem.hisat2_index
    File hisat2_index = BuildHisat2.hisat2_index
    File hisat2_snp_haplotype_splicing_index = BuildHisat2SnpHaplotypeSplicing.hisat2_index
  }
}
