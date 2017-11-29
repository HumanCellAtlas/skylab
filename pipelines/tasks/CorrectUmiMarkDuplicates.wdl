
task CorrectUmiMarkDuplicates {
  File bam_input

  command {
    java -Xmx3g -jar /usr/picard/picard.jar SortSam \
      SORT_ORDER=coordinate \
      I=${bam_input} \
      O=sorted.bam

    java -Xmx3g -jar /usr/picard/picard.jar UmiAwareMarkDuplicatesWithMateCigar \
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
    docker: "humancellatlas/picard:2.10.10"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 100 HDD"
  }
  
  output {
    File bam_output = "duplicates_marked.bam"
    File umi_metrics = "umi_metrics.txt"
    File duplicate_metrics = "duplicate_metrics.txt"
  }
} 