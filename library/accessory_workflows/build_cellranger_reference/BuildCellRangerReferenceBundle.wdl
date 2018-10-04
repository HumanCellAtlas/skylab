version 1.0


workflow BuildCellRangerRef {
  meta {
    description: "Accessory workflow for builing the reference for CellRanger pipeline."
  }

  input {
    String ref_fasta
    String gtf_file
    String ref_name
    String ref_source

    String docker = "quay.io/humancellatlas/secondary-analysis-cellranger:v1.0.0"
    String memory = "30 GB"
    Int boot_disk_size_gb = 12
    Int disk_space = 50
    Int cpu = 8 
  }

  parameter_meta {
    ref_fasta: "Genome reference in fasta format"
    gtf_file: "Gene annotation file in gtf format"
    disk_space: "(optional) the amount of disk space (GB) to provision for this task"
    ref_name: "The expected file name of the output reference"
    docker: "(optional) the docker image containing the runtime environment for this task"
    boot_disk_size_gb: "(optional) the amount of space (GB) of boot disk to provision for this task. Note: due to the large boot disk space the cellranger image requires, this number should be at least greater than 10."
    cpu: "(optional) the number of cpus to provision for this task"
    memory: "(optional) the amount of memory (MB) to provision for this task"
  }

  call BuildCellRangerReference {
    input:
      ref_fasta = ref_fasta,
      gtf_file = gtf_file,
      ref_name = ref_name,
      docker = docker,
      memory = memory,
      boot_disk_size_gb = boot_disk_size_gb,
      disk_space = disk_space,
      cpu = cpu,
      ref_source = ref_source
  }

  output {
    File cellranger_ref = BuildCellRangerReference.cellRangerRef
  }
}

task BuildCellRangerReference {
  input {
    String ref_fasta
    String gtf_file
    String ref_name
    String ref_source
    String docker
    String memory
    Int boot_disk_size_gb
    Int disk_space
    Int cpu
  }

  command <<<
    # This block will do the following works:
    # - Download the gz files of fasta and gtf data according to the URL in the workflow JSON input
    # - Gunzip them and get the name of the files without ".gz"
    # - Filter the GTF file with cellranger
    # - Make the reference with cellranger

    set -euo pipefail

    ref_fasta_filename="$(basename "~{ref_fasta}")"
    ref_fasta_filename="${ref_fasta_filename%.gz}"
    wget -q -O - "~{ref_fasta}"| zcat > "$ref_fasta_filename"

    wget ~{gtf_file}
    gtf_filename="$(basename "~{gtf_file}")"
    gtf_filename="${gtf_filename%.gz}"
    wget -q -O - "~{gtf_file}"| zcat > "$gtf_filename"

    # The reference building process follows the instructions from
    # https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
    if [[ "~{ref_source}"=="Ensembl" ]]; then
      cellranger mkgtf "$gtf_filename" "filtered.$gtf_filename" \
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
    elif [[ "~{ref_source}"=="Gencode" ]]; then
      cellranger mkgtf "$gtf_filename" "filtered.$gtf_filename" \
                    --attribute=gene_type:protein_coding \
                    --attribute=gene_type:lincRNA \
                    --attribute=gene_type:antisense_RNA \
                    --attribute=gene_type:IG_LV_gene \
                    --attribute=gene_type:IG_V_gene \
                    --attribute=gene_type:IG_V_pseudogene \
                    --attribute=gene_type:IG_D_gene \
                    --attribute=gene_type:IG_J_gene \
                    --attribute=gene_type:IG_J_pseudogene \
                    --attribute=gene_type:IG_C_gene \
                    --attribute=gene_type:IG_C_pseudogene \
                    --attribute=gene_type:TR_V_gene \
                    --attribute=gene_type:TR_V_pseudogene \
                    --attribute=gene_type:TR_D_gene \
                    --attribute=gene_type:TR_J_gene \
                    --attribute=gene_type:TR_J_pseudogene \
                    --attribute=gene_type:TR_C_gene
    else
      exit 1;
    fi
    cellranger mkref --genome=GRCh38 \
                     --fasta="$ref_fasta_filename" \
                     --genes="filtered.$gtf_filename" \
                     --nthreads=$(getconf _NPROCESSORS_ONLN)
    tar cvf "~{ref_name}.tar" GRCh38/

  >>>

  runtime {
    bootDiskSizeGb: boot_disk_size_gb
    disks: "local-disk ${disk_space} HDD"
    docker: docker
    cpu: cpu
    memory: memory
  }

  output {
    File cellRangerRef = "${ref_name}.tar"
  }
}
