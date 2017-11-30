## paired-end alignment
task HISAT2PE {
  File hisat2_ref
  File fq1
  File fq2
  String ref_name
  String output_name
  String sample_name
  command {
    tar -zxvf "${hisat2_ref}"
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -1 ${fq1} \
      -2 ${fq2} \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_name}.log \
      --met-file ${output_name}.hisat2.met.txt --met 5 \
      -p 4 -S ${output_name}.sam
    samtools sort -@ 4 -O bam -o "${output_name}.bam" "${output_name}.sam" 
  }
  runtime {
    docker:"humancellatlas/hisat2:2-2.1.0"
    memory:"5 GB"
    disks: "local-disk 100 HDD"
    cpu: "4"
  }
  output {
    File logfile = "${output_name}.log"
    File metfile = "${output_name}.hisat2.met.txt"
    File output_bam = "${output_name}.bam"
  }
}
## run hisat2 for rsem
## increase gap alignment penalty to avoid gap alignment
## --mp 1,1 --np 1 --score-min L,0,-0.1 is default paramesters when rsem runs alignment by using bowtie2
## no-splice-alignment otherwise rsem will report partial alignment error.
task HISAT2rsem {
  File hisat2_ref
  File fq1
  File fq2
  String ref_name
  String output_name
  String sample_name
  command {
    tar -zxvf "${hisat2_ref}" 
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -1 ${fq1} \
      -2 ${fq2} \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file ${output_name}.log \
      --met-file ${output_name}.hisat2.met.txt --met 5 \
      -k 10 \
      --secondary \
      --no-mixed \
      --no-softclip \
      --no-discordant \
      --rdg 99999999,99999999 \
      --rfg 99999999,99999999 \
      --no-spliced-alignment \
      -p 4 -S ${output_name}.sam 
    samtools view -bS "${output_name}.sam" >"${output_name}.bam"
  }
  runtime {
    docker:"humancellatlas/hisat2:2-2.1.0"
    memory:"5 GB"
    disks: "local-disk 50 HDD"
    cpu: "4"
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
  command {
    tar -zxvf "${hisat2_ref}"
    hisat2 -t \
      -x ${ref_name}/${ref_name} \
      -U ${fq} \
      --rg-id=${sample_name} --rg SM:${sample_name} --rg LB:${sample_name} \
      --rg PL:ILLUMINA --rg PU:${sample_name} \
      --new-summary --summary-file "${output_name}.log" \
      --met-file ${output_name}.hisat2.met.txt --met 5 \
      -p 4 -S ${output_name}.sam
      samtools sort -@ 4 -O bam -o "${output_name}.bam" "${output_name}.sam"
  }
  runtime {
    docker:"humancellatlas/hisat2:2-2.1.0"
    memory:"5 GB"
    disks: "local-disk 25 HDD"
    cpu: "4"
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

  command {
    tar -zxvf "${hisat2_ref}"
    hisat2-inspect --ss --snp \
       -s ${ref_name}/${ref_name} > hisat2_inspect.log
    }
  runtime {
    docker:"humancellatlas/hisat2:2-2.1.0"
    memory:"3 GB"
    disks: "local-disk 25 HDD"
    cpu: "1"
  }
  output {
    File logfile ="hisat2_inspect.log"
  }
}
