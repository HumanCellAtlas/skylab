task BuildCellRangerReference {
  String ref_fasta
  String gtf_file
  String ref_name

  String docker = "quay.io/humancellatlas/secondary-analysis-cellranger"
  String memory = "416 GB"
  Int boot_disk_size_gb = 12
  String disk_space = "250"
  Int cpu = 64

  command {
    set -e

    wget ${ref_fasta}
    ref_fasta_tar_filename=$(echo ${ref_fasta} | sed 's:.*/::')
    gunzip $ref_fasta_tar_filename

    ref_fasta_filename=$(echo "$ref_fasta_tar_filename" | rev | cut -c4- | rev)

    wget ${gtf_file}
    gtf_tar_filename=$(echo ${gtf_file} | sed 's:.*/::')
    gunzip $gtf_tar_filename

    gtf_filename=$(echo "$gtf_tar_filename" | rev | cut -c4- | rev)

    cellranger mkgtf $gtf_filename "filtered.$gtf_filename" \
                     --attribute=gene_biotype:protein_coding \
                     --attribute=gene_biotype:lincRNA \
                     --attribute=gene_biotype:antisense \
                     --attribute=gene_biotype:IG_LV_gene \
                     --attribute=gene_biotype:IG_V_gene \
                     --attribute=gene_biotype:IG_V_pseudogene \
                     --attribute=gene_biotype:IG_D_gene \
                     --attribute=gene_biotype:IG_J_gene \
                     --attribute=gene_biotype:IG_J_pseudogene \
                     --attribute=gene_biotype:IG_C_gene \
                     --attribute=gene_biotype:IG_C_pseudogene \
                     --attribute=gene_biotype:TR_V_gene \
                     --attribute=gene_biotype:TR_V_pseudogene \
                     --attribute=gene_biotype:TR_D_gene \
                     --attribute=gene_biotype:TR_J_gene \
                     --attribute=gene_biotype:TR_J_pseudogene \
                     --attribute=gene_biotype:TR_C_gene


    cellranger mkref --genome=GRCh38 \
                     --fasta="$ref_fasta_filename" \
                     --genes="filtered.$gtf_filename" \
                     --nthreads=$(getconf _NPROCESSORS_ONLN)

    tar cvf ${ref_name}.tar GRCh38/
  }
  runtime {
    docker: docker
    memory: memory
    bootDiskSizeGb: boot_disk_size_gb
    disks: "local-disk " + disk_space + " HDD"
    cpu: cpu
  }
  output {
    File cellRangerRef = "${ref_name}.tar"
  }
}

workflow BuildCellRangerRef {
  String fasta
  String gtf
  String ref_name

  call BuildCellRangerReference {
    input:
      ref_fasta = fasta,
      gtf_file = gtf,
      ref_name = ref_name
  }
  output {
    File cellranger_ref = BuildCellRangerReference.cellRangerRef
  }
}
