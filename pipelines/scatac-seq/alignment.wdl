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
            input_bam = AlignPairedEnd.aligned_bam,
            output_snap_basename = 'output.snap',
            genome_name='mm10',
            genome_file_size=genome_file_size
    }
    call SnapCellByBin {
    	 input:
		    snap_input=SnapPre.output_snap,
            bin_size_list = "5000 10000"
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
        String docker_image = "hisplan/snaptools:latest"
    }

    Float input_size = size(input_fastq1, "GiB") + size(input_fastq2, 'GiB') + size(input_reference,'GiB')

    command {
    	# TODO: unzip the reference bundle
        set -euo pipefail
        declare -r TEMP_DIR=`mktemp -d tmpdir_XXXXXX`
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
            --tmp-folder=$TEMP_DIR \
            --overwrite=TRUE \
            --if-sort=True
    }

    output {
        File aligned_bam = output_bam
    }

    runtime {
        docker: docker_image
        cpu: num_threads
        memory: "16 GB"
        disks: "local-disk " + ceil(3.5 * (if input_size < 1 then 1 else input_size )) + " HDD"
    }
}

task SnapPre {
    input {
        File input_bam
        String output_snap_basename
        String genome_name
        File genome_size_file
String docker_image = "hisplan/snaptools:latest"
    }
    command {
        set -euo pipefail
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
    runtime {
        docker: docker_image
        cpu: num_threads
        memory: "16 GB"
        disks: "local-disk 150 HDD"
    }
}

task snapCellByBin {
     input {
        File snap_input
        String bin_size_list
        String docker_image = "hisplan/snaptools:latest"
     }
     command {
        set -euo pipefail
        # Check if this work, because we are just mutating the file
        snaptools snap-add-bmat  \
            --snap-file=~{snap_input}  \
            --bin-size-list ~{bin_size_list}  \
            --verbose=True
     }
     output {
        File output_snap = snap_input
     }
    runtime {
        docker: docker_image
        cpu: num_threads
        memory: "16 GB"
        disks: "local-disk 150 HDD"
    }
}
