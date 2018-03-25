task FastqToUBam {
  File fastq_file
  String sample_id
  String fastq_suffix = ""

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500
  Int cpu = select_first([opt_cpu, 1])
  # estimate that bam is approximately equal in size to fastq, add 20% buffer
  Int disk = select_first([opt_disk, ceil(size(fastq_file, "G") * 2.2)])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    fastq_file: "input fastq file"
    sample_id: "name of sample matching this file, inserted into read group header"
    fastq_suffix: "a suffix to add to the fastq file; useful with mangled file IDs, since picard requires that the file end in .gz or it will not detect the gzipping."
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    # Adds fastq_suffix if it is passed
    if [ ! -z "${fastq_suffix}" ];
    then
        mv "${fastq_file}" "${fastq_file}""${fastq_suffix}"
    fi

    java -Xmx${command_mem_mb}m -jar /usr/picard/picard.jar FastqToSam \
      FASTQ="${fastq_file}""${fastq_suffix}" \
      SORT_ORDER=unsorted \
      OUTPUT=bamfile.bam \
      SAMPLE_NAME="${sample_id}"
  }
  
  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }
  
  output {
    File bam_output = "bamfile.bam"
  }
}
