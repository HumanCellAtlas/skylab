task StarAlignBamSingleEnd {
  File bam_input
  File tar_star_reference

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-40ead6e"])
  Int machine_mem_mb = select_first([opt_memory_gb, ceil(size(tar_star_reference, "G")) + 6]) * 1000
  Int cpu = select_first([opt_cpu, 16])
  # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
  Int disk = select_first([opt_disk, ceil((size(tar_star_reference, "G") * 2) + (size(bam_input, "G") * 2.2))])
  Int preemptible = select_first([opt_preemptible, 0])

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_input: "unaligned bam file containing genomic sequence, tagged with barcode information"
    tar_star_reference: "star reference tarball"
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    # prepare reference
    mkdir genome_reference
    tar -xf "${tar_star_reference}" -C genome_reference --strip-components 1
    rm "${tar_star_reference}"

    STAR \
      --runMode alignReads \
      --runThreadN ${cpu} \
      --genomeDir genome_reference \
      --readFilesIn "${bam_input}" \
      --outSAMtype BAM Unsorted \
      --outSAMattributes All \
      --outFilterMultimapNmax 1 \
      --outSAMunmapped Within \
      --outSAMprimaryFlag AllBestScore \
      --readFilesType SAM SE \
      --readFilesCommand samtools view -h \
      --runRNGseed 777
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} SSD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File bam_output = "Aligned.out.bam"
    File alignment_log = "Log.final.out"
  }

}
