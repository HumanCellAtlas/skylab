
task Mkref {
  # use kallisto index to generate the kallisto index for use with other functions in this
  # module

  File transcriptome_fasta  # fasta file containing transcripts ONLY (NOT a full genome!)
  Int k  # size of kmer in index; must be an odd number.

  command {
    kallisto index \
      --index kallisto.idx \
      -k ${k} \
      "${transcriptome_fasta}"
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-kallisto:v0.2.2-0.43.1"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk 150 HDD"
  }

  output {
    File index = "kallisto.idx"
  }
}


task QuantSingleEnd {
  # run kallisto quant on single-end fastq data from a single cell, defining the fragment length
  # mean and fragment length standard deviation because kallisto cannot infer it from the data
  # without paired reads.
  # kallisto quant pseudo-aligns reads to the transcriptome, and using the identified k-equivalence
  # classes, runs a quantification algorithm to estimate TPM for each transcript.


  File r1  # forward read
  File index  # kallisto index, output of Kallisto.Mkref task
  Int fragment_length_mean = 425  # normal value for 10x distribution
  Int fragment_length_sd = 125  # normal value for 10x distribution

  command {
    kallisto quant \
      --index "${index}" \
      --output-dir . \
      --bootstrap-samples 100 \
      --single \
      --fragment-length ${fragment_length_mean} \
      --sd ${fragment_length_mean} \
      --threads 1 \
      --pseudobam \
      "${r1}" | samtools view -bS - > pseudo.bam

  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-kallisto:v0.2.2-0.43.1"
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    File abundance_h5 = "abundance.h5"
    File abundance_tsv = "abundance.tsv"
    File log = "run_info.json"
    File bam = "pseudo.bam"
  }
}


task QuantPairedEnd {
  # run kallisto quant on paired fastq reads from a single cell
  # kallisto quant pseudo-aligns to the transcriptome, and using the identified k-equivalence
  # classes, runs a quantification algorithm to estimate TPM for each transcript.

  File r1  # forward read
  File index  # kallisto index, output of Kallisto.Mkref task
  File r2  # reverse read

  command {
    kallisto quant \
      --index "${index}" \
      --output-dir . \
      --bootstrap-samples 100 \
      --threads 1 \
      --pseudobam \
      "${r1}" \
      "${r2}" | samtools view -bS - > pseudo.bam
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-kallisto:v0.2.2-0.43.1"
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    File abundance_h5 = "abundance.h5"
    File abundance_tsv = "abundance.tsv"
    File log = "run_info.json"
    File bam = "pseudo.bam"
  }
}

task QuantPairedEndNoBam {
  # run kallisto quant on paired fastq reads from a single cell
  # kallisto quant pseudo-aligns to the transcriptome, and using the identified k-equivalence
  # classes, runs a quantification algorithm to estimate TPM for each transcript.

  File r1  # forward read
  File r2
  File index  # kallisto index, output of Kallisto.Mkref task
  String sample_name
  
  command {
    kallisto quant \
      --index "${index}" \
      --output-dir . \
      --bootstrap-samples 100 \
      --threads 4 \
      "${r1}" "${r2}" 
    mv abundance.h5 "${sample_name}.abundance.h5"
    mv abundance.tsv "${sample_name}.abundance.tsv"
    mv run_info.json "${sample_name}.run_info.json"
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-kallisto:v0.2.2-0.43.1"
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    File abundance_h5 = "${sample_name}.abundance.h5"
    File abundance_tsv = "${sample_name}.abundance.tsv"
    File log = "${sample_name}.run_info.json"
  }
}


task PseudoSingleEnd {
  # run the kallisto pseudo-alignment algorithm on single-ended fastq data from a single cell to
  # generate the number of reads that map to each k-equivalence class or "transcript compatibility
  # class".

  File r1  # forward read
  File index  # kallisto index, output of Kallisto.Mkref task
  Int fragment_length_mean = 425  # normal value for 10x distribution
  Int fragment_length_sd = 125  # normal value for 10x distribution

  command {
    kallisto pseudo \
      --index "${index}" \
      --output-dir . \
      --single \
      --fragment-length ${fragment_length_mean} \
      --sd ${fragment_length_sd} \
      --threads 1 \
      --pseudobam \
      "${r1}" | samtools view -bS - > pseudo.bam
  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-kallisto:v0.2.2-0.43.1"
    cpu: 4  # note that only 1 thread is supported by pseudobam
    memory: "16 GB"
    disks: "local-disk 100 HDD"
  }

  output {
    File pseudo_ec = "pseudoalignments.ec"
    File pseudo_tsv = "pseudoalignments.tsv"
    File run_log = "run_info.json"
    File bam = "pseudo.bam"
  }
}


task PseudoSingleEndUMI {
  # run the kallisto pseudo-alignment algorithm on single-ended fastq data from a single cell to
  # generate the number of unique UMIs that map to each k-equivalence class or "transcript
  # compatibility class".

  File umi  # file containing extracted umis, one per line (not fastq format)
  File r1  # forward read
  File index  # kallisto index, output of Kallisto.Mkref task

  # estimate disk requirements
  Int disks = ceil((size(r1, "G") + size(umi, "G")) * 2.2 + size(index, "G"))

  command {
    echo "cell1 ${umi} ${r1}" >> batch.txt
    kallisto pseudo \
      --index "${index}" \
      --output-dir . \
      --single \
      --umi \
      --threads 4 \
      --batch batch.txt

  }

  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-kallisto:v0.2.2-0.43.1"
    cpu: 4
    memory: "16 GB"
    disks: "local-disk ${disks} HDD"
  }

  output {
    File pseudo_cells = "matrix.cells"
    File pseudo_ec = "matrix.ec"
    File pseudo_tsv = "matrix.tsv"
    File run_log = "run_info.json"
  }
}
