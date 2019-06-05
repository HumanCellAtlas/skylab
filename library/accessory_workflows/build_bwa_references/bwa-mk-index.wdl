version 1.0

workflow BuildBWARef {
  input {
      String ref_name
      File reference_fasta
      File chrom_sizes_file
  }

  call BuildBWAreference {
     input:
	    ref_name = ref_name,
            reference_fasta = reference_fasta,
            chrom_sizes_file = chrom_sizes_file
	}

	output {
	    File referenceBundle = BuildBWAreference.referenceBundle
	}
}

task BuildBWAreference {
     input {
        String ref_name ## name of the tar.bz2 files without the suffix
        File reference_fasta
        File chrom_sizes_file
     }

     command <<<
        mkdir genome
        mv ~{chrom_sizes_file} genome/chrom.sizes
        mv ~{reference_fasta} genome/genome.fa
        bwa index genome/genome.fa
        tar cvf - genome/ > ~{reference_fasta}.tar
     >>>

     runtime {
     	 docker: "hisplan/snaptools"
	     memory: "96GB"
	     disks: "local-disk 100 HDD"
	     cpu: "4"
     }

     output {
     	    File referenceBundle = "~{reference_fasta}.tar"
     }
}
