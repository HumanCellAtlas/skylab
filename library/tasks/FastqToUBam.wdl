task FastqToUBam {
  File fastq_file
  String sample_id
  String fastq_suffix = ""

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-picard:v0.2.2-2.10.10"
  Int machine_mem_mb = 3500
  # give the command 500MB of overhead
  Int command_mem_mb = machine_mem_mb - 500
  Int cpu = 1
  # estimate that bam is approximately equal in size to fastq, add 20% buffer
  Int disk = ceil(size(fastq_file, "G") * 2.2)
  Int preemptible = 0

  meta {
    description: "Converts a fastq file into an unaligned bam file."
  }

  parameter_meta {
    fastq_file: "input fastq file"
    sample_id: "name of sample matching this file, inserted into read group header"
    fastq_suffix: "a suffix to add to the fastq file; useful with mangled file IDs, since picard requires that the file end in .gz or it will not detect the gzipping."
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
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
