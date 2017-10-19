
task StarMkref {
  File fasta_file  # genome annotation
  File annotation_file  # gtf annotation file

  # optionally, use a different version of star
  String star_docker_image = "humancellatlas/star:2.5.3a-40ead6e"
  
  command {
    mkdir genome
    STAR \
      --runMode genomeGenerate \
      --runThreadN $(nproc) \
      --genomeDir genome \
      --genomeFastaFiles "${fasta_file}" \
      --sjdbGTFfile "${annotation_file}"

    echo "tarring genome" && tar -cf genome.tar genome/*
  }
  
  runtime {
    docker: "${star_docker_image}"
    cpu: 16
    memory: "30 GB"
    disks: "local-disk 100 SSD"
  }
  
  output {
    File genome = "genome.tar"
  }
}
