task SortAndCorrectUmiMarkDuplicates {
  File bam_input

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"])
  Int machine_mem_mb = select_first([opt_memory_gb, 7]) * 1000
  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500
  Int cpu = select_first([opt_cpu, 1])
  # mark duplicates swaps a large amount of data to disk, has high disk requirements.
  Int disk = select_first([opt_disk, ceil(size(bam_input, "G") * 6)])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_input: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar SortSam \
      SORT_ORDER=coordinate \
      I=${bam_input} \
      O=sorted.bam

    # recover disk space
    rm ${bam_input}

    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar UmiAwareMarkDuplicatesWithMateCigar \
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
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "duplicates_marked.bam"
    File umi_metrics = "umi_metrics.txt"
    File duplicate_metrics = "duplicate_metrics.txt"
  }
} 
