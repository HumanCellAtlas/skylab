task MergeSortBamFiles {
  Array[File] bam_inputs
  String sort_order

  Int compression_level = 5

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"])
  Int machine_mem_mb = select_first([opt_memory_gb, 7]) * 1000
  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500
  Int cpu = select_first([opt_cpu, 1])
  # default to 500GB of space
  Int disk = select_first([opt_disk, 500])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_inputs: "Merges Sam/Bam files"
    sort_order: "sort order of output bam"
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    java -Dsamjdk.compression_level=${compression_level} -Xms${command_mem_mb}m -jar /usr/gitc/picard.jar \
      MergeSamFiles \
      USE_THREADING=true \
      SORT_ORDER=${sort_order} \
      INPUT=${sep=' INPUT=' bam_inputs} \
      OUTPUT=merged.bam \
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  output {
    File output_bam = "merged.bam"
  }
}
