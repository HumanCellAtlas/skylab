task RSEMExpression {
  File trans_aligned_bam
  File rsem_genome
  String rsem_out

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  Int cpu = select_first([opt_cpu, 4])
  # use provided disk number or dynamically size on our own, with 20GB of additional disk
  Int disk = select_first([opt_disk, ceil(size(trans_aligned_bam, "GB") + size(rsem_genome, "GB") + 20)])
  Int preemptible = select_first([opt_preemptible, 5])

  meta {
    description: "JISHU HELP AGAINNNNNNNNNNNNN"
  }

  parameter_meta {
    trans_aligned_bam: ""
    rsem_genome: ""
    opt_docker: ""
    opt_memory_gb: ""
    opt_cpu: ""
    opt_disk: ""
    opt_preemptible: ""
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
      "${rsem_out}" 
  }

  runtime {
    docker: docker
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File rsem_gene = "${rsem_out}.genes.results"
    File rsem_isoform = "${rsem_out}.isoforms.results"
    File rsem_time = "${rsem_out}.time"
    File rsem_cnt = "${rsem_out}.stat/${rsem_out}.cnt"
    File rsem_model = "${rsem_out}.stat/${rsem_out}.model"
    File rsem_theta = "${rsem_out}.stat/${rsem_out}.theta"
  }
}

