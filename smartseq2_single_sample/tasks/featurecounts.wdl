task FeatureCountsUniqueMapping {
  File aligned_bam
  File gtf
  String fc_out
  
  command {
    featureCounts -T 4  -s 0 -t exon -g gene_id -p -B -C -a "${gtf}" -o "${fc_out}.gene.unq.counts.txt" "${aligned_bam}"
    featureCounts -T 4  -s 0 -t exon -g transcript_id -p -B -C -a "${gtf}" -o "${fc_out}.transcript.unq.counts.txt" "${aligned_bam}"
    featureCounts -T 4  -s 0 -t exon -g exon_id -p -B -C -a "${gtf}" -o "${fc_out}.exon.unq.counts.txt" "${aligned_bam}"
  }
  runtime {
    docker:"humancellatlas/star_dev:2.5.3a"
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
    cpu: "1"
  }
  output {
    File genes = "${fc_out}.gene.unq.counts.txt"
    File exons = "${fc_out}.exon.unq.counts.txt"
    File trans = "${fc_out}.transcript.unq.counts.txt"
  }
}

task FeatureCountsMultiMapping {
  File aligned_bam
  File gtf
  String fc_out

  command {
    featureCounts -T 4 -s 0 -t exon -g gene_id -p -M -O -a "${gtf}" -o "${fc_out}.gene.mult.counts.txt" "${aligned_bam}"
    featureCounts -T 4 -s 0 -t exon -g transcript_id -p -M -O -a "${gtf}" -o "${fc_out}.transcript.mult.counts.txt" "${aligned_bam}"
    featureCounts -T 4 -s 0 -t exon -g exon_id -p -M  -O -a "${gtf}" -o "${fc_out}.exon.mult.counts.txt" "${aligned_bam}"
  }
  runtime {
    docker: "humancellatlas/star_dev:2.5.3a"
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
    cpu: "1"
  }
  output {
    File genes = "${fc_out}.gene.mult.counts.txt"
    File exons = "${fc_out}.exon.mult.counts.txt"
    File trans = "${fc_out}.transcript.mult.counts.txt"
  }

}

