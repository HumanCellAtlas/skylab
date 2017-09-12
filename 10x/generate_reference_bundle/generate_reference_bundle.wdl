
task CellRangerMkref {
  File fasta
  File gtf
  String reference_name

  # todo do we need to fix this attribute for our gtf? (doubt it)
  command {
    cellranger mkgtf ${gtf} filtered.gtf --attribute=gene_biotype:protein_coding
    cellranger mkref --genome=${reference_name} --fasta=${fasta} --genes=filtered.gtf
    tar -cf ${reference_name}.tar ${reference_name}/*
  }

  runtime {
    docker: "marcusczi/cellranger_clean:latest"
    memory: "30 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    File reference_bundle = "${reference_name}.tar"
  }

}


workflow GenerateReferenceBundle {
  File fasta
  File gtf
  String reference_name

  call CellRangerMkref {
    input:
      fasta = fasta,
      gtf = gtf,
      reference_name = reference_name
  }

  output {
    File reference_bundle = CellRangerMkref.reference_bundle
  }

}