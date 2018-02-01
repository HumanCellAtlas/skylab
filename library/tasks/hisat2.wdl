## paired-end alignment
## run HISAT2 to genome reference with dedault parameters
## --seed to fix pseudo-random number and in order to produce deterministics results
## -k --secondary to output multiple mapping reads. --keep 10 will output up to 10 multiple mapping reads, which is default in HISAT2
task HISAT2PE {
  File hisat2_ref
  File fq1  # gz forward fastq file
  File fq2  # gz reverse fastq file
  String ref_name
  String output_name
  String sample_name
  Int disk_size
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
      -p 4 -S ${output_name}.sam 
    samtools sort -@ 4 -O bam -o "${output_name}.bam" "${output_name}.sam" 
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory:"5 GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "4"
    preemptible: 5
  }
  output {
    File logfile = "${output_name}.log"
    File metfile = "${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
}
## run hisat2 for rsem
## increase gap alignment penalty to avoid gap alignment
## --mp 1,1 --np 1 --score-min L,0,-0.1 is default paramesters when rsem runs alignment by using bowtie2/Bowtie
## --mp 1,1 and --np 1 will reduce mismatching penalty to 1 for all.
## with no-splice-alignment no-softclip no-mixed options on, HISAT2 will only output concordant alignment without soft-cliping
## --rdg 99999999,99999999 and --rfg 99999999,99999999 will give an infinity penalty to alignment with indel.As results
## no indel/gaps in alignments
 
task HISAT2rsem {
  File hisat2_ref
  File fq1
  File fq2
  String ref_name
  String output_name
  String sample_name
  Int disk_size
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
      -p 4 -S ${output_name}.sam 
    samtools view -bS  "${output_name}.sam" > "${output_name}.bam"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory:"5 GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "4"
    preemptible: 5
  }
  output {
    File logfile = "${output_name}.log"
    File metfile = "${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
}
## single end alignment
task HISAT2SE {
  File hisat2_ref
  File fq
  String ref_name
  String output_name
  String sample_name
  Int disk_size
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
      -p 4 -S ${output_name}.sam
    samtools sort -@ 4 -O bam -o "${output_name}.bam" "${output_name}.sam"
  }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory:"5 GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "4"
    preemptible: 5
  }
  output {
    File logfile ="${output_name}.log"
    File metfile ="${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
 }
## run inspection on hisat2 built reference index
task hisat2_inspect_index {
  File hisat2_ref
  String ref_name
  Int disk_size
  command {
    set -e
    tar -zxvf "${hisat2_ref}"
    hisat2-inspect --ss --snp \
       -s ${ref_name}/${ref_name} > hisat2_inspect.log
    }
  runtime {
    docker:"quay.io/humancellatlas/secondary-analysis-hisat2:v0.2.2-2-2.1.0"
    memory:"3 GB"
    disks: "local-disk " + disk_size + " HDD"
    cpu: "1"
    preemptible: 5
  }
  output {
    File logfile ="hisat2_inspect.log"
  }
}
