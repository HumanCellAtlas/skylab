task HISAT2PairedEnd {
  File hisat2_ref
  File fq1
  File fq2
  String ref_name
  String output_name
  String sample_name

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"])
  Int machine_mem_mb = select_first([opt_memory_gb, 5]) * 1000
  Int cpu = select_first([opt_cpu, 4])
  # use provided disk number or dynamically size on our own, 10 is our zipped fastq -> bam conversion with 50GB of additional disk
  Int disk = select_first([opt_disk, ceil((size(fq1, "GB") + size(fq2, "GB") * 10) + size(hisat2_ref, "GB") + 50)])
  Int preemptible = select_first([opt_preemptible, 5])

  meta {
    description: "JISHUUUUUUUUUUUUUUUU"
  }

  parameter_meta {
    hisat2_ref: ""
    fq1: "gz forward fastq file"
    fq2: "gz reverse fastq file"
    ref_name: ""
    output_name: ""
    sample_name: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    # Note that files MUST be gzipped or the module will not function properly
    # This will be addressed in the future either by a change in how Hisat2 functions or a more
    # robust test for compression type.

    set -e

    # fix names if necessary.
    if [ "${fq1}" != *.fastq.gz ]; then
        FQ1=${fq1}.fastq.gz
        mv ${fq1} ${fq1}.fastq.gz
    else
        FQ1=${fq1}
    fi
    if [ "${fq2}" != *.fastq.gz ]; then
        FQ2=${fq2}.fastq.gz
        mv ${fq2} ${fq2}.fastq.gz
    else
        FQ2=${fq2}
    fi

   tar -zxvf "${hisat2_ref}"
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -1 $FQ1 \
      -2 $FQ2 \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_name}.log \
      --met-file ${output_name}.hisat2.met.txt --met 5 \
      --seed 12345 \
      -k 10 \
      --secondary \
      -p ${cpu} -S ${output_name}.sam
    samtools sort -@ ${cpu} -O bam -o "${output_name}.bam" "${output_name}.sam"
  }

  runtime {
    docker: docker
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file = "${output_name}.log"
    File met_file = "${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
}

task HISAT2RSEM {
  File hisat2_ref
  File fq1
  File fq2
  String ref_name
  String output_name
  String sample_name

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"])
  Int machine_mem_mb = select_first([opt_memory_gb, 5]) * 1000
  Int cpu = select_first([opt_cpu, 4])
  # use provided disk number or dynamically size on our own, 10 is our zipped fastq -> bam conversion with 50GB of additional disk
  Int disk = select_first([opt_disk, ceil((size(fq1, "GB") + size(fq2, "GB") * 10) + size(hisat2_ref, "GB") + 50)])
  Int preemptible = select_first([opt_preemptible, 5])

  meta {
    description: "JISHUUUUUUUUUUUUUUUU"
  }

  parameter_meta {
    hisat2_ref: ""
    fq1: "gz forward fastq file"
    fq2: "gz reverse fastq file"
    ref_name: ""
    output_name: ""
    sample_name: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e

    # fix names if necessary.
    if [ "${fq1}" != *.fastq.gz ]; then
        FQ1=${fq1}.fastq.gz
        mv ${fq1} ${fq1}.fastq.gz
    else
        FQ1=${fq1}
    fi

    if [ "${fq2}" != *.fastq.gz ]; then
        FQ2=${fq2}.fastq.gz
        mv ${fq2} ${fq2}.fastq.gz
    else
        FQ2=${fq2}
    fi

    tar -zxvf "${hisat2_ref}"
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -1 $FQ1 \
      -2 $FQ2 \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_name}.log \
      --met-file ${output_name}.hisat2.met.txt --met 5 \
      -k 10 \
      --mp 1,1 \
      --np 1 \
      --score-min L,0,-0.1 \
      --secondary \
      --no-mixed \
      --no-softclip \
      --no-discordant \
      --rdg 99999999,99999999 \
      --rfg 99999999,99999999 \
      --no-spliced-alignment \
      --seed 12345 \
      -p ${cpu} -S ${output_name}.sam
    samtools view -bS  "${output_name}.sam" > "${output_name}.bam"
  }

  runtime {
    docker: docker
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file = "${output_name}.log"
    File met_file = "${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
}

task HISAT2SingleEnd {
  File hisat2_ref
  File fq
  String ref_name
  String output_name
  String sample_name

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"])
  Int machine_mem_mb = select_first([opt_memory_gb, 5]) * 1000
  Int cpu = select_first([opt_cpu, 4])
  # use provided disk number or dynamically size on our own, 10 is our zipped fastq -> bam conversion with 50GB of additional disk
  Int disk = select_first([opt_disk, ceil((size(fq, "GB") * 10) + size(hisat2_ref, "GB") + 50)])
  Int preemptible = select_first([opt_preemptible, 5])

  meta {
    description: "JISHUUUUUUUUUUUUUUUU"
  }

  parameter_meta {
    hisat2_ref: ""
    fq: ""
    ref_name: ""
    output_name: ""
    sample_name: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e
    tar -zxvf "${hisat2_ref}"
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -U ${fq} \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file "${output_name}.log" \
      --met-file ${output_name}.hisat2.met.txt --met 5 \
      --seed 12345 \
      -p ${cpu} -S ${output_name}.sam
    samtools sort -@ ${cpu} -O bam -o "${output_name}.bam" "${output_name}.sam"
  }

  runtime {
    docker: docker
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file ="${output_name}.log"
    File met_file ="${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
}

task HISAT2InspectIndex {
  File hisat2_ref
  String ref_name

  # runtime optional arguments
  String? opt_docker
  Int? opt_memory_gb
  Int? opt_cpu
  Int? opt_disk
  Int? opt_preemptible

  # runtime values
  String docker = select_first([opt_docker, "quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"])
  Int machine_mem_mb = select_first([opt_memory_gb, 3]) * 1000
  Int cpu = select_first([opt_cpu, 1])
  # use provided disk number or dynamically size on our own, with 10GB of additional disk
  Int disk = select_first([opt_disk, ceil(size(hisat2_ref, "GB") + 10)])
  Int preemptible = select_first([opt_preemptible, 5])

  meta {
    description: "JISHUUUUUUUUUUUUUUUU"
  }

  parameter_meta {
    hisat2_ref: ""
    ref_name: ""
    opt_docker: "optionally provide a docker to run in"
    opt_memory_gb: "optionally provide how much memory to provision"
    opt_cpu: "optionally provide how many cpus to provision"
    opt_disk: "optionally provide how much disk to provision"
    opt_preemptible: "optionally provide how many preemptible attempts"
  }

  command {
    set -e
    tar -zxvf "${hisat2_ref}"
    hisat2-inspect --ss --snp \
       -s ${ref_name}/${ref_name} > hisat2_inspect.log
  }

  runtime {
    docker: docker
    memory: machine_mem_mb + " MB"
    disks: "local-disk " + disk + " HDD"
    cpu: cpu
    preemptible: preemptible
  }

  output {
    File log_file ="hisat2_inspect.log"
  }
}
