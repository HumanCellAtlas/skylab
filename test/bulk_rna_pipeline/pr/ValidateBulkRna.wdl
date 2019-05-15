task ValidateBulkRna {
      File anno_flagstat
      File annobam
      File genome_flagstat
      File genomebam
      File quants
      File mad_qc_metrics
      File rna_qc_metrics
      File rsem_genes_detected
      String expected_anno_flagstat_hash
      String expected_annobam_hash
      String expected_genome_flagstat_hash
      String expected_genomebam_hash
      String expected_quants_hash
      String expected_mad_qc_metrics_hash
      String expected_rna_qc_metrics_hash
      String expected_rsem_genes_detected_hash

  command <<<

    # catch intermittent failures
    set -eo pipefail

    # calculate hash, but ignore run-specific bam headers
    annobam_hash=$(samtools view "${annobam}" | md5sum | awk '{print $1}')
    genomebam_hash=$(samtools view "${genomebam}" | md5sum | awk '{print $1}')

    anno_flagstat_hash=$(cat "${anno_flagstat}" | md5sum | awk '{print $1}')
    genome_flagstat_hash=$(cat "${genome_flagstat}" | md5sum | awk '{print $1}')
    quants_hash=$(cat "${quants}" | md5sum | awk '{print $1}')
    mad_qc_metrics_hash=$(cat "${mad_qc_metrics}" | md5sum | awk '{print $1}')
    rna_qc_metrics_hash=$(cat "${rna_qc_metrics}" | md5sum | awk '{print $1}')
    rsem_genes_detected_hash=$(cat "${rsem_genes_detected}" | md5sum | awk '{print $1}')

    if [ "$annobam_hash" != "${expected_annobam_hash}" ]; then
      >&2 echo "annobam_hash ($annobam_hash) did not match expected hash (${expected_annobam_hash})"
      fail=true
    fi

    if [ "$genomebam_hash" != "${expected_genomebam_hash}" ]; then
      >&2 echo "genomebam_hash ($genomebam_hash) did not match expected hash (${expected_genomebam_hash})"
      fail=true
    fi

    if [ "$anno_flagstat_hash" != "${expected_anno_flagstat_hash}" ]; then
      >&2 echo "anno_flagstat_hash ($anno_flagstat_hash) did not match expected hash (${expected_anno_flagstat_hash})"
      fail=true
    fi

    if [ "$genome_flagstat_hash" != "${expected_genome_flagstat_hash}" ]; then
      >&2 echo "genome_flagstat_hash ($genome_flagstat_hash) did not match expected hash (${expected_genome_flagstat_hash})"
      fail=true
    fi

    if [ "$quants_hash" != "${expected_quants_hash}" ]; then
      >&2 echo "quants_hash ($quants_hash) did not match expected hash (${expected_quants_hash})"
      fail=true
    fi

    if [ "$mad_qc_metrics_hash" != "${expected_mad_qc_metrics_hash}" ]; then
      >&2 echo "mad_qc_metrics_hash ($mad_qc_metrics_hash) did not match expected hash (${expected_mad_qc_metrics_hash})"
      fail=true
    fi

    if [ "$rna_qc_metrics_hash" != "${expected_rna_qc_metrics_hash}" ]; then
      >&2 echo "rna_qc_metrics_hash ($rna_qc_metrics_hash) did not match expected hash (${expected_rna_qc_metrics_hash})"
      fail=true
    fi

    if [ "$rsem_genes_detected_hash" != "${expected_rsem_genes_detected_hash}" ]; then
      >&2 echo "rsem_genes_detected_hash ($rsem_genes_detected_hash) did not match expected hash (${expected_rsem_genes_detected_hash})"
      fail=true
    fi

    if [ $fail == "true" ]; then exit 1; fi

  >>>

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
  }
}
