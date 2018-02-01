
task CorrectUmiMarkDuplicates {
  File bam_input
  # mark duplicates swaps a large amount of data to disk, has high disk requirements.
  Int estimated_required_disk = ceil(size(bam_input, "G") * 6)

  command {
    java -Xmx7g -jar /usr/picard/picard.jar SortSam \
      SORT_ORDER=coordinate \
      I=${bam_input} \
      O=sorted.bam

    # recover disk space
    rm ${bam_input}

    java -Xmx7g -jar /usr/picard/picard.jar UmiAwareMarkDuplicatesWithMateCigar \
      MAX_EDIT_DISTANCE_TO_JOIN=1 \
      UMI_METRICS_FILE=umi_metrics.txt \
      UMI_TAG_NAME=UR \
      ASSIGNED_UMI_TAG=UB \
      BARCODE_TAG=CB \
      METRICS_FILE=duplicate_metrics.txt \
      OUTPUT=duplicates_marked.bam \
      INPUT=sorted.bam
  }
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
    cpu: 1
    memory: "7.5 GB"
    disks: "local-disk ${estimated_required_disk} HDD"
  }
  
  output {
    File bam_output = "duplicates_marked.bam"
    File umi_metrics = "umi_metrics.txt"
    File duplicate_metrics = "duplicate_metrics.txt"
  }
} 