task BuildBWAreference {
  String ref_name ## name of the tar.bz2 files without the suffix
     File reference_fasta

     command {
        mkdir genome
        mv ${reference_fasta} genome/genome.fa
        bwa index genome/genome.fa
        tar cvf - genome/ > ${reference_fasta}.tar
     }

     runtime {
     	 docker: "hisplan/snaptools"
	     memory: "96GB"
	     disks: "local-disk 100 HDD"
	     cpu: "4"
     }

     output {
     	    File referenceBundle = "${reference_fasta}.tar"
     }
}

workflow BuildBWARef {
  String ref_name
  File reference_fasta

  call BuildBWAreference {
     input:
	    ref_name = ref_name,
		reference_fasta = reference_fasta
	}

	output {
	       File referenceBundle = BuildBWAreference.referenceBundle
	}
}