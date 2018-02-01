## RSEM will estimate gene/isoform quantification
## --bam: input is bam file
## -p 4: run on multiple threads
## --time: report running time
## --seed: report deterministic results
## --calc-pme will do posterior mean estimation
## --single-cell-prior in pme estimation, use Dirichlet(0.1) as the prior 
## which encourage the sparsity of the expression levels
## of note, the --single-cell-prior has not been tested yet and there are
## more investigation or improvement requested.  
task RsemExpression {
  File trans_aligned_bam
  File rsem_genome
  String rsem_out
  Int disk_size
  command {
    set -e
  
    tar -xvf ${rsem_genome}
    rsem-calculate-expression \
      --bam \
      --paired-end \
       -p 4 \
      --time --seed 555 \
      --calc-pme \
      --single-cell-prior \
      ${trans_aligned_bam} \
      rsem/rsem_trans_index  \
      "${rsem_out}" 
  }
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-rsem:v0.2.2-1.3.0"
    memory: "3.75 GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "4"
    preemptible: 5
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

