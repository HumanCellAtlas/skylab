task RSEMExpression {
  File trans_aligned_bam
  File rsem_genome
  String output_basename

  # runtime values
  String docker = "quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0"
  Int machine_mem_mb = 3500
  Int cpu = 4
  # use provided disk number or dynamically size on our own, with 20GB of additional disk
  Int disk = ceil(size(trans_aligned_bam, "GB") + size(rsem_genome, "GB") + 20)
  Int preemptible = 5

  meta {
    description: "This task will quantify gene expression matrix by using RSEM. The output include gene-level and isoform-level results."
  }

  parameter_meta {
    trans_aligned_bam: ""
    rsem_genome: ""
    output_basename: "basename used for output files"
    docker: "optionally provide a docker image"
    machine_mem_mb: "optionally provide how much memory(MB) to provision"
    cpu: "optionally provide how many cpus to provision"
    disk: "optionally provide how much disk to provision"
    preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e
  
    tar -xvf ${rsem_genome}
    rsem-calculate-expression \
      --bam \
      --paired-end \
       -p ${cpu} \
      --time --seed 555 \
      --calc-pme \
      --single-cell-prior \
      ${trans_aligned_bam} \
      rsem/rsem_trans_index  \
      "${output_basename}"
  }

  runtime {
    docker: docker
    memory: "${machine_mem_mb} MB"
    disks: "local-disk ${disk} HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File rsem_gene = "${output_basename}.genes.results"
    File rsem_isoform = "${output_basename}.isoforms.results"
    File rsem_time = "${output_basename}.time"
    File rsem_cnt = "${output_basename}.stat/${output_basename}.cnt"
    File rsem_model = "${output_basename}.stat/${output_basename}.model"
    File rsem_theta = "${output_basename}.stat/${output_basename}.theta"
  }
}

