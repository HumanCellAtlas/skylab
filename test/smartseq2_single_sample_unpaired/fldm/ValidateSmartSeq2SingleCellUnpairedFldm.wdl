task ValidateSmartSeq2SingleCellUnpairedFldm {
      File alignment_summary_metrics
      File gene_counts
      String expected_gene_counts_hash

      String dollar='$'

  command <<<

    # catch intermittent failures
    set -xeo pipefail

    # Verify that the PCT_PF_ALIGNED_READS is reasonable.
    # Strip out comments and blank lines and get the genotype concordance
    PCT_PF_READS_ALIGNED=$(grep -e CATEGORY -e UNPAIRED  "${alignment_summary_metrics}" |
      sed '/#/d' |
      sed '/^\s*$/d' |
      awk -v col=PCT_PF_READS_ALIGNED \
      'NR==1{for(i=1;i<=NF;i++){if($i==col){c=i;break}}} NR==2{print $c}')


    # TODO - make this a reasonable value NOT 0.01
    if awk 'BEGIN{exit ARGV[1]>=ARGV[2]}' "$PCT_PF_READS_ALIGNED" "0.01"
    then
      # Less than threshold, fail.
      echo "% PF Aligned Reads below threshold. Failure"
      exit 1
    fi

    # calculate hashes; awk is used to extract the hash from the md5sum output that contains both
    # a hash and the filename that was passed. We parse the first 7 columns because a bug in RSEM
    # makes later columns non-deterministic.
    gene_counts_hash=$(cut -f 1-7 "${gene_counts}" | md5sum | awk '{print $1}')

    if [ "$gene_counts_hash" != "${expected_gene_counts_hash}" ]; then
      >&2 echo "counts_hash ($gene_counts_hash) did not match expected hash (${expected_gene_counts_hash})"
      fail=true
    fi

    if [ $fail == "true" ]; then exit 1; fi

  >>>
  
  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk 10 HDD"
  }
}
