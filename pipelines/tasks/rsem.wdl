task RsemExpression {
  File trans_aligned_bam
  File rsem_genome
  String rsem_out
  
  command {
    tar -xvf ${rsem_genome}
    echo "Aligning fastqs and calculating expression"
    rsem-calculate-expression --bam --paired-end ${trans_aligned_bam} rsem/rsem_trans_index  "${rsem_out}" -p 4 --ci-memory 2000 --time 
    ## parse gene expected_count out
    cut -f 1,4,5 "${rsem_out}.genes.results" >"${rsem_out}.gene.expected_counts"
  }
  runtime {
    docker: "humancellatlas/rsem:1.3.0"
    memory: "3.75 GB"
    disks: "local-disk 50 HDD"
    cpu: "4"
  }
  output {
    File rsem_gene = "${rsem_out}.genes.results"
    File rsem_isoform = "${rsem_out}.isoforms.results"
    File rsem_gene_count = "${rsem_out}.gene.expected_counts"
    File rsem_time = "${rsem_out}.time"
    File rsem_cnt = "${rsem_out}.stat/${rsem_out}.cnt"
    File rsem_model = "${rsem_out}.stat/${rsem_out}.model"
    File rsem_theta = "${rsem_out}.stat/${rsem_out}.theta"
  }
}

