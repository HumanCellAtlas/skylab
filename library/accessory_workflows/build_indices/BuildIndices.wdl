version 1.0

struct References {
  File genome_fa
  File annotation_gtf
}

task GetReferences {
  input {
    String gtf_version
    String organism
    String organism_prefix
  }

  String ftp_path = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_~{organism}/release_~{gtf_version}"
  String genome_fa = "GRC~{organism_prefix}38.primary_assembly.genome.fa"
  String annotation_gtf = "gencode.v~{gtf_version}.primary_assembly.annotation.gtf"

  command <<<
    set -eo pipefail

    ## download fasta
    wget ~{ftp_path}/~{genome_fa}.gz
    gunzip ~{genome_fa}.gz

    ## download gtf file
    wget ~{ftp_path}/~{annotation_gtf}.gz
    gunzip ~{annotation_gtf}.gz
  >>>

  output {
      References references = object {
      genome_fa: genome_fa,
      annotation_gtf: annotation_gtf
    }
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-1.0.0"
    disks: "local-disk 10 HDD"
  }
}

task BuildStar {
  input {
    String gtf_version
    String organism
    References references
  }

  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String star_index_name = "~{ref_name}.tar"


  command <<<
    set -eo pipefail

    mkdir star
    STAR --runMode genomeGenerate \
      --genomeDir star \
      --genomeFastaFiles ~{references.genome_fa} \
      --sjdbGTFfile ~{references.annotation_gtf} \
      --sjdbOverhang 100 \
      --runThreadN 16

    tar -cvf ~{star_index_name} star
  >>>

  output {
    File star_index = star_index_name
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-1.0.0"
    memory: "50 GiB"
    disks :"local-disk 100 HDD"
    cpu:"16"
  }
}

task BuildRsem {
  input {
    String gtf_version
    String organism
    References references
  }

  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String rsem_index_name = "~{ref_name}.tar"

  command {
    set -eo pipefail
    mkdir rsem
    rsem-prepare-reference --gtf ~{references.annotation_gtf} --bowtie ~{references.genome_fa} rsem/rsem_trans_index
    tar -cvf ~{rsem_index_name} rsem/
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0"
    memory: "10 GB"
    disks: "local-disk 100 HDD"
  }
  output {
    File rsem_index = rsem_index_name
  }
}

task BuildHisat2FromRsem {
  input {
    File rsem_index
  }

  String rsem_reference_name = basename(rsem_index, ".tar")
  String ref_name = "hisat2_from_rsem_~{rsem_reference_name}"
  String hisat2_index_name = "~{ref_name}.tar.gz"

  command {

    # extract rsem index
    tar -xf ~{rsem_index}

    # build index
    hisat2-build -p 8 rsem/rsem_trans_index.idx.fa ~{ref_name}
    mkdir ~{ref_name}
    mv ./*.ht2 ~{ref_name}
    tar -zcvf ~{hisat2_index_name} ~{ref_name}
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "8 GiB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }

  output {
    File hisat2_index = hisat2_index_name
  }
}

task BuildHisat2 {
  # This method builds a HISAT2 index directly from the downloaded GTF and Fasta files on GENCODE
  # It does not include dbsnp, whereas the below method does.

  input {
    String gtf_version
    String organism
    File genome_fa
  }

  String ref_name = "hisat2_primary_gencode_~{organism}_v~{gtf_version}"
  String hisat2_index_name = "~{ref_name}.tar.gz"

  command {
    set -eo pipefail

    # build index
    hisat2-build -p 8 ~{genome_fa} ~{ref_name}
    mkdir ~{ref_name}
    mv ./*.ht2 ~{ref_name}
    tar -zcvf ~{hisat2_index_name} ~{ref_name}
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory: "8 GiB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }

  output {
    File hisat2_index = hisat2_index_name
  }
}

task BuildHisat2SnpHaplotypeSplicing {
  # This version includes SNP, haplotype, and
  input {
    String organism
    String genome_short_string
    References references

    String gtf_version
    String dbsnp_version
  }

  String ref_name = "star_primary_gencode_~{organism}_v~{gtf_version}"
  String snp_file = "snp~{dbsnp_version}Common.txt"
  String hisat2_index_name = "~{ref_name}.tar.gz"

  command <<<

    HISAT2_DIR=/opt/tools/hisat2-2.1.0

    # Compressed fasta required here
    gzip ~{references.genome_fa}

    # download snp file
    wget http://hgdownload.cse.ucsc.edu/goldenPath/~{genome_short_string}/database/~{snp_file}.gz
    gunzip ~{snp_file}.gz

    # extract snps, splice sites, and exon information
    $HISAT2_DIR/hisat2_extract_snps_UCSC.py ~{references.genome_fa}.gz ~{snp_file} genome
    $HISAT2_DIR/hisat2_extract_splice_sites.py ~{references.annotation_gtf} > genome.ss
    $HISAT2_DIR/hisat2_extract_exons.py ~{references.annotation_gtf} > genome.exon

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
    tar -zcvf ~{hisat2_index_name} ~{ref_name}

  >>>

  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.3.0-2-2.1.0"
    memory: "240 GiB"
    disks: "local-disk 100 HDD"
    cpu: "16"
  }
  output {
    File hisat2_index = hisat2_index_name
  }
}

task BuildPicardRefFlat {
  input {
    References references
  }

  String refflat = basename(references.annotation_gtf, ".gtf") + ".refflat.txt"

  command {
    set -eo pipefail

    gtfToGenePred -genePredExt -geneNameAsName2  ~{references.annotation_gtf} refflat.tmp.txt

    paste <(cut -f 12 refflat.tmp.txt) <(cut -f 1-10 refflat.tmp.txt) > ~{refflat}

  }

  runtime {
    docker: "quay.io/humancellatlas/gtf_to_genepred:v0.0.0"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }

  output {
      File refflat = refflat
  }
}

task BuildIntervalList {
  input {
    References references
  }

  String interval_list = basename(references.annotation_gtf, ".gtf") + ".interval_list"

  command <<<
    set -eo pipefail


    # index the fasta file

    samtools faidx ~{references.genome_fa}
    cut -f1,2 ~{references.genome_fa}.fai > sizes.genome

    awk -F '\t'  '{  printf "@SQ\tSN:%s\tLN:%s\n", $1, $2 }' sizes.genome  >> ~{interval_list}

    grep 'gene_type "rRNA"' ~{references.annotation_gtf} |
        awk '$3 == "transcript"' |
    cut -f1,4,5,7,9 | \
    perl -lane '
        /transcript_id "([^"]+)"/ or die "no transcript_id on $.";
        print join "\t", (@F[0,1,2,3], $1)
    ' | \
    sort -k1V -k2n -k3n  >> ~{interval_list}

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-umitools:0.0.1"
    memory: "8 GB"
    disks: "local-disk 100 HDD"
    cpu: "8"
  }

  output {
      File interval_list = "${interval_list}"
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

  parameter_meta {
    gtf_version: "the actual number of gencode, ex.  27"
    organism: "Either 'human' or 'mouse'"
    organism_prefix: "Either 'h' or 'm'"
    genome_short_string: "e.g. hg38, mm10"
    dbsnp_version: "integer num, ex 150"
  }

  call GetReferences {
    input:
      gtf_version = gtf_version,
      organism = organism,
      organism_prefix = organism_prefix
  }

  call BuildPicardRefFlat   {
     input:
        references = GetReferences.references
  }

  call BuildIntervalList  {
      input:
        references = GetReferences.references
  }

  call BuildStar {
    input:
      gtf_version = gtf_version,
      organism = organism,
      references = GetReferences.references
  }

  call BuildRsem {
    input:
      gtf_version = gtf_version,
      organism = organism,
      references = GetReferences.references
  }

  call BuildHisat2FromRsem {
    input:
      rsem_index = BuildRsem.rsem_index
  }

  call BuildHisat2 {
    input:
      gtf_version = gtf_version,
      organism = organism,
      genome_fa = GetReferences.references.genome_fa
  }

  call BuildHisat2SnpHaplotypeSplicing {
    input:
      gtf_version = gtf_version,
      organism = organism,
      genome_short_string = genome_short_string,
      dbsnp_version = dbsnp_version,
      references = GetReferences.references
  }

  output {
    File star_index = BuildStar.star_index
    File rsem_index = BuildRsem.rsem_index
    File hisat2_from_rsem_index = BuildHisat2FromRsem.hisat2_index
    File hisat2_index = BuildHisat2.hisat2_index
    File hisat2_snp_haplotype_splicing_index = BuildHisat2SnpHaplotypeSplicing.hisat2_index
  }
}
