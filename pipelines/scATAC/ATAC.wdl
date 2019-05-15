version 1.0

workflow ATAC {
 input {
   #fastq inputs
   File fastq_gzipped_input_read1
   File fastq_gzipped_input_read2

   # trimming options
   Int min_length
   Int quality_cutoff
   String adapter_seq_read1
   String adapter_seq_read2

   # BWA: pre built reference, default read group name and # of cores
   File tar_bwa_reference
   String read_group_id = "foo"
   String read_group_sample_name = "bar"
   Int bwa_cpu = 16

   # genome name and genome size file
   String genome_name
   File genome_size_file

   # filtering options
   Int min_map_quailty
   Int max_fragment_length

   # output prefix/base name for all intermediate files and pipeline outputs
   String output_base_name
 }

  call TrimAdapters {
    input:
      fastq_input_read1 = fastq_gzipped_input_read1,
      fastq_input_read2 = fastq_gzipped_input_read2,
      min_length = min_length,
      quality_cutoff = quality_cutoff,
      adapter_seq_read1 = adapter_seq_read1,
      adapter_seq_read2 = adapter_seq_read2,
      output_base_name = output_base_name
  }

  call BWAPairedEndALignment {
    input:
      fastq_input_read1 = TrimAdapters.fastq_trimmed_adapter_output_read1,
      fastq_input_read2 = TrimAdapters.fastq_trimmed_adapter_output_read2,
      tar_bwa_reference = tar_bwa_reference,
      read_group_id = read_group_id,
      read_group_sample_name = read_group_sample_name,
      cpu = bwa_cpu,
      output_base_name = output_base_name
  }

  call SamToBam {
    input:
      sam_input = BWAPairedEndALignment.sam_aligned_output,
      output_base_name = output_base_name
  }

  call MakeCompliantBAM as MakeCompliantAlignedBAM {
    input:
      bam_input = SamToBam.bam_output,
      output_base_name = output_base_name + ".aligned"
  }

  call Sort as SortCoordinateOrder {
    input:
      bam_input = SamToBam.bam_output,
      sort_order = "coordinate",
      output_base_name = output_base_name
  }

  call FilterMarkDuplicates {
    input:
      bam_input = SortCoordinateOrder.bam_sort_output,
      output_base_name = output_base_name
  }

  call FilterMinMapQuality{
    input:
      bam_input = FilterMarkDuplicates.bam_remove_dup_output,
      min_map_quality = min_map_quailty,
      output_base_name = output_base_name
  }

  call FilterMaxFragmentLength {
    input:
      bam_input = FilterMinMapQuality.bam_filter_mapq_output,
      max_fragment_length = max_fragment_length,
      output_base_name = output_base_name
  }

  call FilterMitochondrialReads {
    input:
      bam_input = FilterMaxFragmentLength.bam_filter_fragment_length_output,
      output_base_name = output_base_name
  }

  call MakeCompliantBAM as MakeCompliantChrMBAM {
     input:
      bam_input = FilterMitochondrialReads.bam_chrM_reads_output,
      output_base_name = output_base_name + ".chrM_reads"
  }

  call Sort as SortQueryName {
    input:
      bam_input = FilterMitochondrialReads.bam_no_chrM_reads_output,
      sort_order = "queryname",
      output_base_name = output_base_name
  }

  call MakeCompliantBAM as MakeCompliantFilteredAndSortedBAM {
    input:
      bam_input = SortQueryName.bam_sort_output,
      output_base_name = output_base_name + ".filtered_and_sorted"
  }

  call SnapPre {
    input:
      bam_input= SortQueryName.bam_sort_output,
      output_base_name = output_base_name,
      genome_name = genome_name,
      max_fragment_length = max_fragment_length,
      genome_size_file = genome_size_file
  }

  call SnapCellByBin {
    input:
      snap_input=SnapPre.snap_file_output,
      bin_size_list = "5000 10000"
  }

  output {
    File bam_aligned_compliant_output = MakeCompliantAlignedBAM.compliant_bam_output
    File bam_chrM_reads_compliant_output = MakeCompliantChrMBAM.compliant_bam_output
    File bam_filtered_and_sorted_compliant_output = MakeCompliantFilteredAndSortedBAM.compliant_bam_output
    File snap_qc_output = SnapPre.snap_qc_output
    File snap_output = SnapCellByBin.snap_output
  }
}

task TrimAdapters {
  input {
    Int min_length
    Int quality_cutoff
    File fastq_input_read1
    File fastq_input_read2
    String adapter_seq_read1
    String adapter_seq_read2
    String output_base_name
    String docker_image = "quay.io/broadinstitute/cutadapt:1.18"
  }

  # runtime requirements based upon input file size
  Float input_size = size(fastq_input_read1, "GiB") + size(fastq_input_read2, "GiB")
  Int disk_size = ceil(2 * (if input_size < 1 then 1 else input_size))

  # output names for trimmed reads
  String fastq_trimmed_adapter_output_name_read1 = output_base_name + ".R1.trimmed_adapters.fastq.gz"
  String fastq_trimmed_adapter_output_name_read2 = output_base_name + ".R2.trimmed_adapters.fastq.gz"

  # using cutadapt to trim off sequence adapters
  command {
    set -euo pipefail

    # fastq's, "-f", -A for paired adapters read 2"
    cutadapt \
      -f fastq \
      --minimum-length ~{min_length} \
      --quality-cutoff ~{quality_cutoff} \
      --adapter ~{adapter_seq_read1} \
      -A ~{adapter_seq_read2} \
      --output ~{fastq_trimmed_adapter_output_name_read1} \
      --paired-output ~{fastq_trimmed_adapter_output_name_read2} \
      ~{fastq_input_read1} ~{fastq_input_read2}
  }

  # use docker image for given tool cutadapat
  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: 1
    memory: "3.75 GiB"
  }

  output {
    File fastq_trimmed_adapter_output_read1 = fastq_trimmed_adapter_output_name_read1
    File fastq_trimmed_adapter_output_read2 = fastq_trimmed_adapter_output_name_read2
    File monitoring_log = "monitoring.log"
  }
}

# TODO:
# update docker image to be tool specific
task BWAPairedEndALignment {
  input {
    File fastq_input_read1
    File fastq_input_read2
    File tar_bwa_reference
    String read_group_id
    String read_group_sample_name
    Int cpu
    String output_base_name
    String docker_image = "quay.io/humancellatlas/snaptools:0.0.1"
  }

  # runtime requirements based upon input file size
  Float input_size = size(fastq_input_read1, "GiB") + size(fastq_input_read2, "GiB") + size(tar_bwa_reference, "GiB")
  Int disk_size = ceil(3.25 * (if input_size < 1 then 1 else input_size))

  String sam_aligned_output_name = output_base_name + ".aligned.sam"

  # sort with samtools
  command {
    set -euo pipefail

    # prepare reference
    declare -r REF_DIR=$(mktemp -d genome_referenceXXXXXX)
    tar -xf "~{tar_bwa_reference}" -C $REF_DIR --strip-components 1
    rm "~{tar_bwa_reference}"

    # align w/ BWA: -t for number of cores
    bwa \
      mem \
      -R "@RG\tID:~{read_group_id}\tSM:~{read_group_sample_name}" \
      -t ~{cpu} \
      $REF_DIR/genome.fa \
      <(zcat ~{fastq_input_read1}) <(zcat ~{fastq_input_read2}) \
      > ~{sam_aligned_output_name}
  }

  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: cpu
    memory: "3.75 GiB"
  }

  output {
    File sam_aligned_output = sam_aligned_output_name
    File monitoring_log = "monitoring.log"
  }
}

# filter bam with a minimum mapping quality
task SamToBam {
  input {
    File sam_input
    String output_base_name
    String docker_image = "quay.io/broadinstitute/samtools:1.9"
  }

  # output name for filtered read
  String bam_output_name = output_base_name + ".bam"

  # runtime requirements based upon input file size
  Int disk_size = ceil(2 * (if size(sam_input, "GiB") < 1 then 1 else size(sam_input, "GiB")))

  command {
    set -euo pipefail

    # converst sam to bam
    samtools view \
      -bhS \
      ~{sam_input} \
      -o ~{bam_output_name}
  }

  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: 1
    memory: "3.75 GiB"
  }

  output {
    File bam_output = bam_output_name
    File monitoring_log = "monitoring.log"
  }
}

task Sort {
  input {
    File bam_input
    String sort_order
    String output_base_name
    String docker_image = "quay.io/broadinstitute/picard:2.18.23"
  }

  # output name for sorted bam
  String bam_sort_output_name = output_base_name + ".sorted." + sort_order + ".bam"

  # runtime requirements based upon input file size
  Int disk_size = ceil(3.25 * (if size(bam_input, "GiB") < 1 then 1 else size(bam_input, "GiB")))

  # sort with samtools
  command {
    set -euo pipefail

    java -jar /picard-tools/picard.jar SortSam \
      INPUT=~{bam_input} \
      SORT_ORDER=~{sort_order} \
      MAX_RECORDS_IN_RAM=300000 \
      OUTPUT=~{bam_sort_output_name}
  }

  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: 1
    memory: "3.75 GiB"
  }

  output {
    File bam_sort_output = bam_sort_output_name
    File monitoring_log = "monitoring.log"
  }
}

# filter bam by removing duplicates
task FilterMarkDuplicates {
  input {
    File bam_input
    String output_base_name
    String docker_image = "quay.io/broadinstitute/picard:2.18.23"
  }

  # output namefor mark duplicates
  String bam_remove_dup_output_name = output_base_name + ".filtered.duplicates.bam"
  String metric_remove_dup_output_name = output_base_name + ".filtered.duplicate_metrics"

  # runtime requirements based upon input file size
  Int disk_size = ceil(2 * (if size(bam_input, "GiB") < 1 then 1 else size(bam_input, "GiB")))

  command {
    set -euo pipefail

    java -jar /picard-tools/picard.jar MarkDuplicates \
      INPUT=~{bam_input} \
      OUTPUT=~{bam_remove_dup_output_name} \
      METRICS_FILE=~{metric_remove_dup_output_name}
  }

  runtime {
     docker: docker_image
     disks: "local-disk " + disk_size + " HDD"
     cpu: 1
     memory: "3.75 GiB"
  }

  output {
    File bam_remove_dup_output = bam_remove_dup_output_name
    File metric_remove_dup_output = metric_remove_dup_output_name
    File monitoring_log = "monitoring.log"
  }
}

# filter bam with a minimum mapping quality
task FilterMinMapQuality {
  input {
    File bam_input
    Int min_map_quality
    String output_base_name
    String docker_image = "quay.io/broadinstitute/samtools:1.9"
  }

  # output name for filtered read
  String bam_filter_mapq_output_name = output_base_name + ".filtered.min_map_quality.bam"

  # runtime requirements based upon input file size
  Int disk_size = ceil(2 * (if size(bam_input, "GiB") < 1 then 1 else size(bam_input, "GiB")))

  command {
    set -euo pipefail

    # filter for a map quality
    # -b output is bam, -h include header, -q reads with mapping quality >=
    samtools view \
      -bh \
      -q~{min_map_quality} \
      ~{bam_input} \
      > ~{bam_filter_mapq_output_name}
  }

  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: 1
    memory: "3.75 GiB"
  }

  output {
    File bam_filter_mapq_output = bam_filter_mapq_output_name
    File monitoring_log = "monitoring.log"
  }
}

# filter bam with a max fragment length
task FilterMaxFragmentLength {
  input {
    File bam_input
    Int max_fragment_length
    String output_base_name
    String docker_image = "broadinstitute/gatk:4.1.2.0"
  }

  # output name for filtered read
  String bam_filter_fragment_length_output_name = output_base_name + ".filtered.max_fragment_length.bam"

  # runtime requirements based upon input file size
  Int disk_size = ceil(2 * (if size(bam_input, "GiB") < 1 then 1 else size(bam_input, "GiB")))

  command {
    set -euo pipefail

    gatk PrintReads \
      --input=~{bam_input} \
      --read-filter FragmentLengthReadFilter --max-fragment-length ~{max_fragment_length} \
      --output=~{bam_filter_fragment_length_output_name}
  }

  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: 1
    memory: "3.75 GiB"
  }

  output {
    File bam_filter_fragment_length_output = bam_filter_fragment_length_output_name
    File monitoring_log = "monitoring.log"
  }
}

# split the bam into read containing only and not containing any mitochondrial reads
task FilterMitochondrialReads {
  input {
    File bam_input
    String output_base_name
    String docker_image = "quay.io/broadinstitute/samtools:1.9"
  }

  # output name for sorted bam
  String bam_chrM_reads_output_name = output_base_name + ".filtered.mitochondrial_reads.bam"
  String bam_no_chrM_reads_output_name = output_base_name +".filtered.no_mitochondrial_reads.bam"

  # runtime requirements based upon input file size
  Int disk_size = ceil(2 * (if size(bam_input, "GiB") < 1 then 1 else size(bam_input, "GiB")))

  # ChrM: mitochondrial chromosome
  command {
    set -euo pipefail

    # create bam index to filter by chromosome
    samtools index -b ~{bam_input}

    # get bam w/o chrM
    list_chrs=`samtools view -H ~{bam_input} \
      | grep chr \
      | cut -f2 \
      | sed 's/SN://g' \
      | grep -v 'chrM\|_'`

    samtools view \
      -bh \
      -f 0x2 \
      ~{bam_input} \
      `echo $list_chrs` \
      -o ~{bam_no_chrM_reads_output_name}

    #get bam with only chrM
    samtools view  \
      -bh \
      ~{bam_input} \
      chrM \
      -o ~{bam_chrM_reads_output_name}
  }

  runtime {
    docker: docker_image
    disks: "local-disk " + disk_size + " HDD"
    cpu: 1
    memory: "3.75 GiB"
  }

  output {
    File bam_no_chrM_reads_output = bam_no_chrM_reads_output_name
    File bam_chrM_reads_output = bam_chrM_reads_output_name
    File monitoring_log = "monitoring.log"
  }
}

task SnapPre {
  input {
    File bam_input
    String output_base_name
    String genome_name
    Int max_fragment_length
    File genome_size_file
    String docker_image = "quay.io/humancellatlas/snaptools:0.0.1"
  }

  String snap_file_output_name = output_base_name + ".snap"
  String snap_qc_output_name = snap_file_output_name + ".qc"


  # TODO:
  # update disk size to be dynamilcally updated based off of input file size

  command {
    set -euo pipefail

    # Does the main counting
    snaptools snap-pre \
      --input-file=~{bam_input} \
      --output-snap=~{snap_file_output_name} \
      --genome-name=~{genome_name} \
      --genome-size=~{genome_size_file} \
      --min-mapq=0  \
      --min-flen=0  \
      --max-flen=~{max_fragment_length}  \
      --keep-chrm=TRUE  \
      --keep-single=TRUE  \
      --keep-secondary=False  \
      --overwrite=True  \
      --max-num=1000000  \
      --min-cov=100  \
      --verbose=True
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "16 GiB"
    disks: "local-disk 150 HDD"
  }

  output {
    File snap_file_output = snap_file_output_name
    File snap_qc_output = snap_qc_output_name
  }
}

task SnapCellByBin {
  input {
    File snap_input
    String bin_size_list
    String docker_image = "quay.io/humancellatlas/snaptools:0.0.1"
  }

  command {
    set -euo pipefail

    # This is mutating the file in-place
    snaptools snap-add-bmat  \
      --snap-file=~{snap_input}  \
      --bin-size-list ~{bin_size_list}  \
      --verbose=True
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "16 GiB"
    disks: "local-disk 150 HDD"
  }

  output {
    File snap_output = snap_input
  }
}

# TODO:
# update docker image to be tool specific
task MakeCompliantBAM {
  input {
    File bam_input
    String output_base_name
    String docker_image = "quay.io/humancellatlas/snaptools:0.0.1"
  }

  Int disk_size = ceil(2.5 * (if size(bam_input, "GiB") < 1 then 1 else size(bam_input, "GiB")))

  String compliant_bam_output_name = output_base_name + ".compliant.bam"

  command {
    makeCompliantBAM.py \
      --input-bam ~{bam_input} \
      --output-bam ~{compliant_bam_output_name}
  }

  runtime {
    docker: docker_image
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File compliant_bam_output = compliant_bam_output_name
  }
}
