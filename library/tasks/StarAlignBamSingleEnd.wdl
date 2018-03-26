task StarAlignBamSingleEnd {
  File bam_input
  File tar_star_reference

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-star:v0.2.2-2.5.3a-40ead6e"
  Int machine_mem_mb = (ceil(size(tar_star_reference, "G")) + 6) * 1000
  Int cpu = 16
  # multiply input size by 2.2 to account for output bam file + 20% overhead, add size of reference.
  Int disk = ceil((size(tar_star_reference, "G") * 2) + (size(bam_input, "G") * 2.2))
  Int preemptible = 0

  meta {
    description: "AMBROSE HALP!!"
  }

  parameter_meta {
    bam_input: "unaligned bam file containing genomic sequence, tagged with barcode information"
    tar_star_reference: "star reference tarball"
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
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
