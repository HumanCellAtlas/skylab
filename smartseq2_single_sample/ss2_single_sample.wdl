task Star {
  File input_fastq_read1
  File input_fastq_read2
  File gtf
  File star_genome

  command {
    tar -xvf ${star_genome}
    STAR  --readFilesIn ${input_fastq_read1} ${input_fastq_read2} \
      --genomeDir ./star \
      --quantMode TranscriptomeSAM \
      --outSAMstrandField intronMotif \
      --genomeLoad NoSharedMemory \
      --sjdbGTFfile ${gtf} \
      --readFilesCommand "zcat" \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate  \
      --outSAMunmapped Within \
      --limitBAMsortRAM 30000000000 
  }
  output {
    File junction_table = "SJ.out.tab"
    File output_bam = "Aligned.sortedByCoord.out.bam"
    File output_bam_trans="Aligned.toTranscriptome.out.bam"
  }
  runtime {
    docker:"humancellatlas/star:2.5.3a"
    memory: "40 GB"
    disks :"local-disk 100 HDD"
  }
}

task FeatureCountsUniqueMapping {
  File aligned_bam
  File gtf
  String fc_out
  
  command {
    featureCounts -s 0 -t exon -g gene_id -p -B -C -a ${gtf} -o "${fc_out}.gene.counts.txt" ${aligned_bam}
    featureCounts -s 0 -t exon -g transcript_id -p -B -C -a ${gtf} -o "${fc_out}.transcripts.counts.txt" ${aligned_bam}
    featureCounts -s 0 -t exon -g exon_id -p -B -C -a ${gtf} -o "${fc_out}.exon.counts.txt" ${aligned_bam}
  }
  runtime {
    docker:"humancellatlas/star:2.5.3a"
    memory: "15 GB"
    disks :"local-disk 50 HDD"
  }
  output {
    File genes="${fc_out}.gene.counts.txt"
    File exons="${fc_out}.exon.counts.txt"
    File trans="${fc_out}.transcripts.counts.txt"
  }
}

task RsemExpression {
  File trans_aligned_bam
  File rsem_genome
  String rsem_out
  
  command {
    tar -xvf ${rsem_genome}
    echo "Aligning fastqs and calculating expression"
    rsem-calculate-expression --bam --paired-end ${trans_aligned_bam} rsem/rsem_trans_index  "${rsem_out}"
    ## parse gene expected_count out
    cut -f 1,4,5 "${rsem_out}.genes.results" >"${rsem_out}.gene.expected_counts"
  }
  runtime {
    docker: "humancellatlas/rsem:v1.3.0"
    memory: "10 GB"
    disks: "local-disk 100 HDD"
  }
  output {
    File rsem_gene="${rsem_out}.genes.results"
    File rsem_transc="${rsem_out}.isoforms.results"
    File rsem_gene_count="${rsem_out}.gene.expected_counts"
   }
}

task CollectAlignmentSummaryMetrics {
  File aligned_bam
  File ref_genome_fasta
  String output_filename
  
  command {
    java -Xmx10g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT=${aligned_bam} \
      OUTPUT="${output_filename}" \
      REFERENCE_SEQUENCE=${ref_genome_fasta} \
      ASSUME_SORTED=true
  }
  output {
    File alignment_metrics ="${output_filename}"
  }
  runtime {
    docker:"humancellatlas/picard:2.10.10"
    memory:"10 GB"
    disks: "local-disk 10 HDD"
  }
}

task CollectRnaSeqMetrics {
  File aligned_bam
  File ref_genome_fasta
  File rrna_interval
  String output_filename
  File ref_flat
  
  command {
    java -Xmx10g -jar /usr/picard/picard.jar  CollectRnaSeqMetrics \
      VALIDATION_STRINGENCY=SILENT \
      REF_FLAT=${ref_flat} \
      RIBOSOMAL_INTERVALS=${rrna_interval} \
      INPUT=${aligned_bam} \
      OUTPUT="${output_filename}" \
      REFERENCE_SEQUENCE=${ref_genome_fasta} \
      ASSUME_SORTED=true \
      STRAND_SPECIFICITY=NONE
  }
  output {
    File rna_metrics="${output_filename}"
  }
  runtime {
    docker:"humancellatlas/picard:2.10.10"
    memory:"10 GB"
    disks: "local-disk 10 HDD"
  }
}

workflow Ss2RsemSingleSample {
  File fastq_read1
  File fastq_read2
  File gtf
  File ref_fasta
  File rrna_interval
  File ref_flat
  String star_genome
  String output_prefix
  String rsem_genome
  
  call Star {
    input:
      input_fastq_read1 = fastq_read1,
      input_fastq_read2 = fastq_read2,
      gtf = gtf,
      star_genome = star_genome
  }
 
  call RsemExpression {
    input:
      trans_aligned_bam=Star.output_bam_trans,
      rsem_genome = rsem_genome,
      rsem_out = output_prefix
  }

  call FeatureCountsUniqueMapping {
    input:
      aligned_bam=Star.output_bam,
      gtf =gtf,
      fc_out=output_prefix
  }

  call CollectRnaSeqMetrics {
    input:
      aligned_bam=Star.output_bam,
      ref_genome_fasta=ref_fasta,
      rrna_interval = rrna_interval,
      output_filename = "${output_prefix}_rna_metrics",
      ref_flat= ref_flat
  }

  call CollectAlignmentSummaryMetrics {
    input:
      aligned_bam = Star.output_bam,
      ref_genome_fasta = ref_fasta,
      output_filename ="${output_prefix}_alignment_metrics"
  }

  output {
    File bam_file=Star.output_bam
    File bam_trans=Star.output_bam_trans
    File rna_metrics=CollectRnaSeqMetrics.rna_metrics
    File aln_metrics=CollectAlignmentSummaryMetrics.alignment_metrics
    File rsem_gene_results=RsemExpression.rsem_gene
    File rsem_isoform_results=RsemExpression.rsem_transc
    File rsem_gene_count=RsemExpression.rsem_gene_count
    File gene_unique_counts=FeatureCountsUniqueMapping.genes
    File exon_unique_counts=FeatureCountsUniqueMapping.exons
    File transcript_unique_counts=FeatureCountsUniqueMapping.trans
  }
}
