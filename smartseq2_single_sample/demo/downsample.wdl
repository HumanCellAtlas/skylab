task DownSampleByChrm {
  File bamfile
  String chrm
  String frac

  command {
    samtools index ${bamfile}
    samtools view -h  -s ${frac}  ${bamfile} ${chrm} \
      |samtools sort -n -O bam \
      |samtools fastq -n - -1 R1.fastq -2 R2.fastq
    gzip R1.fastq R2.fastq
  }
  output {
    File R1 = "R1.fastq.gz"
    File R2 = "R2.fastq.gz"
    }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
    memory: "10 GB"
    disks :"local-disk 10 HDD"
  }
}

workflow DownSampleBam {
  File filename
  String chrN
  String fraction

  call DownSampleByChrm {
    input:
      bamfile = filename,
      chrm = chrN,
      frac =fraction
  }
  output {
    File r1 = DownSampleByChrm.R1
    File r2 = DownSampleByChrm.R2
  }
}
