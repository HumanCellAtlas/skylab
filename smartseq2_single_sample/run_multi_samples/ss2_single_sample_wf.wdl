task Star {
  File input_fastq_read1
  File input_fastq_read2
  File gtf
  File star_genome
  String sample_tag
  String pu_tag
  String id_tag
  String lib_tag


  command {
    tar -xvf ${star_genome}
    STAR  --readFilesIn ${input_fastq_read1} ${input_fastq_read2} \
      --genomeDir ./star \
      --quantMode TranscriptomeSAM \
      --outSAMstrandField intronMotif \
      --runThreadN 8 \
      --genomeLoad NoSharedMemory \
      --sjdbGTFfile ${gtf} \
      --readFilesCommand "zcat" \
      --twopassMode Basic \
      --outSAMtype BAM SortedByCoordinate  \
      --outSAMunmapped Within \
      --limitBAMsortRAM 30000000000 \
      --outSAMattrRGline ID:${id_tag} PL:illumina PU:${pu_tag} SM:${sample_tag} LB:${lib_tag}
  }
  output {
    File junction_table = "SJ.out.tab"
    File output_bam = "Aligned.sortedByCoord.out.bam"
    File output_bam_trans = "Aligned.toTranscriptome.out.bam"
  }
  runtime {
    docker:"humancellatlas/star:2.5.3a"
    memory: "40 GB"
    disks :"local-disk 100 HDD"
    cpu: "8"
  }
}

task FeatureCountsUniqueMapping {
  File aligned_bam
  File gtf
  String fc_out
  
  command {
    featureCounts -s 0 -t exon -g gene_id -p -B -C -a ${gtf} -o "${fc_out}.gene.unq.counts.txt" ${aligned_bam}
    featureCounts  -s 0 -t exon -g transcript_id -p -B -C -a ${gtf} -o "${fc_out}.transcript.unq.counts.txt" ${aligned_bam}
    featureCounts -s 0 -t exon -g exon_id -p -B -C -a ${gtf} -o "${fc_out}.exon.unq.counts.txt" ${aligned_bam}
  }
  runtime {
    docker:"humancellatlas/star:2.5.3a"
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
  output {
    File genes = "${fc_out}.gene.unq.counts.txt"
    File exons = "${fc_out}.exon.unq.counts.txt"
    File trans = "${fc_out}.transcript.unq.counts.txt"
  }
}

task FeatureCountsMultiMapping {
  File aligned_bam
  File gtf
  String fc_out

  command {
    featureCounts  -s 0 -t exon -g gene_id -p -M -O -a ${gtf} -o "${fc_out}.gene.mult.counts.txt" ${aligned_bam}
    featureCounts  -s 0 -t exon -g transcript_id -p -M -O -a ${gtf} -o "${fc_out}.transcript.mult.counts.txt" ${aligned_bam}
    featureCounts  -s 0 -t exon -g exon_id -p -M  -O -a ${gtf} -o "${fc_out}.exon.mult.counts.txt" ${aligned_bam}
  }
  runtime {
    docker: "humancellatlas/star:2.5.3a"
    memory: "4 GB"
    disks: "local-disk 50 HDD"
  }
  output {
    File genes = "${fc_out}.gene.mult.counts.txt"
    File exons = "${fc_out}.exon.mult.counts.txt"
    File trans = "${fc_out}.transcript.mult.counts.txt"
  }
}

task RsemExpression {
  File trans_aligned_bam
  File rsem_genome
  String rsem_out
  
  command {
    tar -xvf ${rsem_genome}
    echo "Aligning fastqs and calculating expression"
    rsem-calculate-expression --bam --paired-end ${trans_aligned_bam} rsem/rsem_trans_index  "${rsem_out}" -p 4 --ci-memory 2000
    ## parse gene expected_count out
    cut -f 1,4,5 "${rsem_out}.genes.results" >"${rsem_out}.gene.expected_counts"
  }
  runtime {
    docker: "humancellatlas/rsem:v1.3.0"
    memory: "4 GB"
    disks: "local-disk 50 HDD"
    cpu: "4"
  }
  output {
    File rsem_gene = "${rsem_out}.genes.results"
    File rsem_transc = "${rsem_out}.isoforms.results"
    File rsem_gene_count = "${rsem_out}.gene.expected_counts"
   }
}

task CollectAlignmentSummaryMetrics {
  File aligned_bam
  File ref_genome_fasta
  String output_filename
  
  command {
    java -Xmx4g -jar /usr/gitc/picard.jar CollectAlignmentSummaryMetrics \
      VALIDATION_STRINGENCY=SILENT \
      METRIC_ACCUMULATION_LEVEL=ALL_READS \
      INPUT=${aligned_bam} \
      OUTPUT="${output_filename}.alignment_metrics" \
      REFERENCE_SEQUENCE=${ref_genome_fasta} \
      ASSUME_SORTED=true
  }
  output {
    File alignment_metrics = "${output_filename}.alignment_metrics"
  }
  runtime {
    docker:"broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
    memory:"4 GB"
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
    java -Xmx4g -jar /usr/gitc/picard.jar  CollectRnaSeqMetrics \
      VALIDATION_STRINGENCY=SILENT \
      REF_FLAT=${ref_flat} \
      RIBOSOMAL_INTERVALS=${rrna_interval} \
      INPUT=${aligned_bam} \
      OUTPUT="${output_filename}.rna_metrics" \
      REFERENCE_SEQUENCE=${ref_genome_fasta} \
      ASSUME_SORTED=true \
      STRAND_SPECIFICITY=NONE
  }
  output {
    File rna_metrics = "${output_filename}.rna_metrics"
  }
  runtime {
    docker:"broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
    memory:"4 GB"
    disks: "local-disk 10 HDD"
  }
}

task CollectDuplicationMetrics {
  File aligned_bam
  String output_filename

  command {
    java -Xmx4g -jar /usr/gitc/picard.jar  MarkDuplicates \
       VALIDATION_STRINGENCY=SILENT  \
       INPUT=${aligned_bam} \
       OUTPUT="${output_filename}.MarkDuplicated.bam" \
       ASSUME_SORTED=true \
       METRICS_FILE="${output_filename}.duplicate_metrics" \
       REMOVE_DUPLICATES=false
  }
  output {
    File dedup_metrics = "${output_filename}.duplicate_metrics"
    File dedup_bamfile = "${output_filename}.MarkDuplicated.bam"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
    memory: "4 GB"
    disks: "local-disk 20 HDD"
  }
}

task CollectInsertMetrics {
  File aligned_bam
  String output_filename

  command {
    java -Xmx4g -jar /usr/gitc/picard.jar CollectInsertSizeMetrics \
      INPUT=${aligned_bam} \
      OUTPUT="${output_filename}.insert_size_metrics" \
      HISTOGRAM_FILE="${output_filename}.insert_size_histogram.pdf" \
  }
  output {
    File insert_metrics = "${output_filename}.insert_size_metrics"
    File histogram = "${output_filename}.insert_size_histogram.pdf"
  }
  runtime {
    docker: "broadinstitute/genomes-in-the-cloud:2.3.1-1504795437"
    memory: "4 GB"
    disks: "local-disk 10 HDD"
  }
}
task CollectMetricsBySample {
  Array[File] input_metrics_list
  String output_filename
  
  command {
    python /tools/scripts/parse_picard.py -I ${write_lines(input_metrics_list)} -T "sample" -O "${output_filename}"
  }
  runtime {
    docker: "humancellatlas/python3-application:0.1.0"
    memory: "2 GB"
    disk: "local-disk 10 HDD"
  }
  output {
    File sample_metrics = "${output_filename}.metrics"
  }
}
workflow Ss2RunSingleSample {
  File fastq_read1
  File fastq_read2
  File gtf
  File ref_fasta
  File rrna_interval
  File ref_flat
  String star_genome
  String rsem_genome
  String output_prefix
  
  call Star {
    input:
      input_fastq_read1 = fastq_read1,
      input_fastq_read2 = fastq_read2,
      gtf = gtf,
      star_genome = star_genome,
      sample_tag = output_prefix,
      pu_tag = output_prefix,
      lib_tag =output_prefix,
      id_tag = output_prefix
  }
 
  call RsemExpression {
    input:
      trans_aligned_bam = Star.output_bam_trans,
      rsem_genome = rsem_genome,
      rsem_out = output_prefix
  }

  call FeatureCountsUniqueMapping {
    input:
      aligned_bam = Star.output_bam,
      gtf = gtf,
      fc_out = output_prefix
  }
  
  call FeatureCountsMultiMapping {
    input:
      aligned_bam=Star.output_bam,
      gtf = gtf,
      fc_out = output_prefix
  }

  
  call CollectRnaSeqMetrics {
    input:
      aligned_bam = Star.output_bam,
      ref_genome_fasta = ref_fasta,
      rrna_interval = rrna_interval,
      output_filename = "${output_prefix}",
      ref_flat= ref_flat
  }

  call CollectAlignmentSummaryMetrics {
    input:
      aligned_bam = Star.output_bam,
      ref_genome_fasta = ref_fasta,
      output_filename = "${output_prefix}"
  }
  
  call CollectDuplicationMetrics {
    input:
      aligned_bam = Star.output_bam,
      output_filename = "${output_prefix}"
  }

  call CollectInsertMetrics {
    input:
      aligned_bam = Star.output_bam,
      output_filename = "${output_prefix}"
    }
  
  call CollectMetricsBySample {
    input:
      input_metrics_list = [CollectAlignmentSummaryMetrics.alignment_metrics,CollectRnaSeqMetrics.rna_metrics,CollectInsertMetrics.insert_metrics,CollectDuplicationMetrics.dedup_metrics],
      output_filename = "${output_prefix}"
  }
  output {
    File bam_file = Star.output_bam
    File bam_trans = Star.output_bam_trans
    File rna_metrics = CollectRnaSeqMetrics.rna_metrics
    File aln_metrics = CollectAlignmentSummaryMetrics.alignment_metrics
    File dedup_metrics = CollectDuplicationMetrics.dedup_metrics
    File dedup_bam = CollectDuplicationMetrics.dedup_bamfile
    File insert_metrics = CollectInsertMetrics.insert_metrics
    File insert_histogram = CollectInsertMetrics.histogram
    File sample_metrics = CollectMetricsBySample.sample_metrics
    File rsem_gene_results = RsemExpression.rsem_gene
    File rsem_isoform_results = RsemExpression.rsem_transc
    File rsem_gene_count = RsemExpression.rsem_gene_count
    File gene_unique_counts = FeatureCountsUniqueMapping.genes
    File exon_unique_counts = FeatureCountsUniqueMapping.exons
    File transcript_unique_counts = FeatureCountsUniqueMapping.trans
    File gene_multi_counts = FeatureCountsMultiMapping.genes
    File exon_multi_counts = FeatureCountsMultiMapping.exons
    File transcript_multi_counts = FeatureCountsMultiMapping.trans
  }
}
