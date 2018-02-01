task htseq_count {
  File aligned_bam
  File gtf 
  String featuretype
  String output_filename
  command {
    htseq-count --format=bam --stranded=no \
      --order=pos \
      --type="${featuretype}" \
      --idattr="${featuretype}_id" \
      "${aligned_bam}" \
      "${gtf}" \
      > "${output_filename}.${featuretype}.htseq.count.txt"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-python3-scientific:0.1.5"
    memory: "7.5 GB"
    disks: "local-disk 100 HDD"
    cpu: "2"
  }
  output {
    File counts = "${output_filename}.${featuretype}.htseq.count.txt"	
  }
}
