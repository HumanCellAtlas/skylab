
task MergeSortBam {
  Array[Array[File]] bam_inputs

  command {
    flat_inputs=$(python -c "print '${sep=' ' bam_inputs}'.replace('[', '').replace(']', '')")
    samtools merge -c -p out.bam $flat_inputs
    samtools sort -o out_sorted.bam out.bam
  }
  
  runtime {
    docker: "humancellatlas/samtools:1.3.1"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 100 HDD"
  }
  
  output {
    File bam_output = "out_sorted.bam"
  }
}
