task MergeSortBamFiles {
  Array[File] bam_inputs
  String sort_order

  Int compression_level = 5

  # runtime values
  String docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
  Int machine_mem_mb = 7500
  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500
  Int cpu = 1
  # default to 500GB of space
  Int disk = 500
  Int preemptible = 0

  meta {
    description: "Merge multiple bam files in the specified sort order"
  }

  parameter_meta {
    bam_inputs: "Merges Sam/Bam files"
    sort_order: "sort order of output bam"
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
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
