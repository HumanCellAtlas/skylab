version 1.0

workflow scATAC {
    input {
        File input_fastq1
        File input_fastq2
        File input_reference
        String output_bam
        String genome_name
        File genome_size_file
    }
    parameter_meta {
        input_fastq1: "read 1 input fastq"
        input_fastq2: "read 2 input fastq"
        input_reference: "tar file with BWA reference"
        output_bam: "output BAM file"
        genome_name: "name of the genome for snap atac"
        genome_size_file: "two column tsv file with contig names and lengths"
    }
    call AlignPairedEnd {
        input:
            input_fastq1 = input_fastq1,
            input_fastq2 = input_fastq2,
            input_reference = input_reference,
            output_bam = output_bam
    }
    call SnapPre {
        input:
            input_bam = AlignPairedEnd.aligned_bam,
            output_snap_basename = 'output.snap',
            genome_name = genome_name,
            genome_size_file = genome_size_file
    }
    call SnapCellByBin {
        input:
            snap_input=SnapPre.output_snap,
            bin_size_list = "5000 10000"
    }
    output {
        File output_snap_qc = SnapPre.output_snap_qc
        File output_snap = SnapCellByBin.output_snap
        File output_aligned_bam = AlignPairedEnd.aligned_bam
    }
}

task AlignPairedEnd {
    input {
        File input_fastq1
        File input_fastq2
        File input_reference
        String reference_unpack_name = "genome/genome.fa"
        String output_bam
        Int min_cov = 0
        String docker_image = "hisplan/snaptools:latest"
    }

    Int num_threads = 1
    Float input_size = size(input_fastq1, "GiB") + size(input_fastq2, "GiB") + size(input_reference, "GiB")

    command {
        set -euo pipefail

        # Make temp directory
        declare -r TEMP_DIR=`mktemp -d tmpdir_XXXXXX`

        # Unpack the reference
        tar xf ~{input_reference}

        # Run snaptools alignment
        snaptools align-paired-end \
            --input-reference=~{reference_unpack_name} \
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

    Int num_threads = 1

    command {
        set -euo pipefail

        # Does the main counting
        snaptools snap-pre \
            --input-file=~{input_bam} \
            --output-snap=~{output_snap_basename} \
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

task SnapCellByBin {
    input {
        File snap_input
        String bin_size_list
        String docker_image = "hisplan/snaptools:latest"
    }

    Int num_threads = 1

    command {
        set -euo pipefail

        # This is mutating the file in-place
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
