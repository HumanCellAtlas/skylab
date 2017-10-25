task Star {
  File input_fastq_read1
  File input_fastq_read2
  File gtf
  File star_genome
  String sample_tag
  String pu_tag
  String id_tag
  String lib_tag

  command {
    tar -xvf ${star_genome}
      
    STAR  --readFilesIn ${input_fastq_read1} ${input_fastq_read2} \
      --genomeDir ./star \
      --quantMode TranscriptomeSAM \
      --outSAMstrandField intronMotif \
      --runThreadN 8 \
      --genomeLoad NoSharedMemory \
      --sjdbGTFfile ${gtf} \
      --readFilesCommand "zcat" \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate  \
      --outSAMunmapped Within \
      --limitBAMsortRAM 30000000000 \
      --outFileNamePrefix "${sample_tag}." \
      --outSAMattrRGline ID:${id_tag} PL:illumina PU:${pu_tag} SM:${sample_tag} LB:${lib_tag}
    
    mkdir logs
    mv ${sample_tag}.Log.* ./logs/
    tar -zcvf  ${sample_tag}.logs.tar.gz logs/
  }
  runtime {
    docker:"humancellatlas/star_dev:2.5.3a"
    memory: "40 GB"
    disks :"local-disk 100 HDD"
    cpu: "8"
  }
  output {
    File junction_table = "${sample_tag}.SJ.out.tab"
    File logs = "${sample_tag}.logs.tar.gz"
    File output_bam = "${sample_tag}.Aligned.sortedByCoord.out.bam"
    File output_bam_trans = "${sample_tag}.Aligned.toTranscriptome.out.bam"
  }
}

