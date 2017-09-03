task BuildStarReference{
  File ref_fasta
  File gtf_file
  
  command {
    mkdir star
    STAR --runMode genomeGenerate \
      --genomeDir star \
      --genomeFastaFiles ${ref_fasta} \
      --sjdbGTFfile ${gtf_file} \
      --sjdbOverhang 100 \
      --runThreadN 16
    tar -cvf star.tar star
  }
  output {
    File starRef = "star.tar"
  }
  runtime {
    docker:"humancellatlas/star_dev:v1"
    memory: "50 GB"
    disks :"local-disk 100 HDD"
    cpu:"16"
  }
}


workflow StarRef {
  File fasta
  File gtf
  
  call BuildStarReference {
    input:
      ref_fasta = fasta,
      gtf_file = gtf
  }
  output {
    File star_ref = BuildStarReference.starRef
  }
}
