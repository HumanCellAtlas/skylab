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
    docker:"humancellatlas/python3-scientific:0.1.1"
    memory: "3.75 GB"
    disks: "local-disk 50 HDD"
    cpu: "1"
   }
  output {
    File counts = "${output_filename}.${featuretype}.htseq.count.txt"	
    }
}
