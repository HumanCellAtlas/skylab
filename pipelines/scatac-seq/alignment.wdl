version 1.0

workflow scATAC {
    input {
        File input_fastq1
        File input_fastq2
        File input_reference
        String output_bam

	String genome_name
	File genome_fize_file

    }
    call AlignPairedEnd {
        input:
            input_fastq1 = input_fastq1,
            input_fastq2 = input_fastq2,
            input_reference = input_reference,
            output_bam = output_bam,
            min_cov = min_cov,
            num_threads = num_threads,
    }
    call SnapPre {
        input:
            input_bam = AlignPairedEnd.output_bam,
	    output_snap_basename = 'output.snap',
	    genome_name='mm10',
	    genome_file_size=genome_file_size
    }
    call SnapCellByBin {
    	 input:
		snap_input=SnapPre.output_snap,
                String bin_size_list = "5000 10000"
    }
}

task AlignPairedEnd {
    input {
        File input_fastq1
        File input_fastq2
        File input_reference
        File output_bam
        Int min_cov=0
        Int num_threads=1
    }

    command {
        set -euo pipefail
	mkdir -p tmp/
        snaptools align-paired-end \
            --input-reference=~{input_reference} \
            --input-fastq1=~{input_fastq1} \
            --input-fastq2=~{input_fastq2} \
            --output-bam=~{output_bam} \
            --aligner=bwa \
            --path-to-aligner=/tools/ \
            --read-fastq-command=zcat \
            --min-cov=~{min_cov} \
            --num-threads=~{num_threads} \
            --tmp-folder=tmp/
            --overwrite=TRUE \
            --if-sort=True
    }

    output {
        File path_output = "~{output_bam}"
    }

    runtime {
        docker: "hisplan/snaptools:latest"
        cpu: 1
        memory: "16 GB"
        disks: "local-disk 150 HDD"
    }

}

task snapPre {
    input {
    	  File input_bam
	  String output_snap_basename
	  String genome_name
	  File genome_size_file
    }
    command {
        set -euo pipefail \
        snaptools snap-pre \
            --input-file=~{input_bam} \
            --output-snap=~{output_snap_basenname} \
            --genome-name=~{genome_name} \
            --genome-size=~{genome_size_file} \
	    --min-mapq=30  \
	    --min-flen=0  \
	    --max-flen=1000  \
	    --keep-chrm=TRUE  \
	    --keep-single=TRUE  \
	    --keep-secondary=False  \
	    --overwrite=True  \
	    --max-num=1000000  \
	    --min-cov=100  \
	    --verbose=True
    }   
    output {
    	   File output_snap = output_snap_basename
	   File output_snap_qc = output_snap_basename + ".qc"
    }    
}

task snapCellByBin {
     input {
     	   File snap_input
	   String bin_size_list
     }
     command {
     	  # Check if this work, because we are just mutating the file
          snaptools snap-add-bmat  \
     	  	    --snap-file=~{snap_input}  \
     		    --bin-size-list ~{bin_size_list}  \
     		    --verbose=True
     }
     output {
     	    File output_snap = snap_input
     }
}
