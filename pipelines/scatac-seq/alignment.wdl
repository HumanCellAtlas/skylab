version 1.0

workflow scATAC {

    input {
        File input_fastq1
        File input_fastq2
        File input_reference
        String output_bam
        String aligner
        String path_aligner
        String read_fastq_command
        Int min_cov
        Int num_threads
        Boolean sort
        String tmp_folder
        Boolean overwrite

    }

    call alignPairedEnd {
        input:
            input_fastq1 = input_fastq1,
            input_fastq2 = input_fastq2,
            input_reference = input_reference,
            output_bam = output_bam,
            aligner = aligner,
            path_aligner = path_aligner,
            read_fastq_command = read_fastq_command,
            min_cov = min_cov,
            num_threads = num_threads,
            sort = sort,
            tmp_folder = tmp_folder,
            overwrite = overwrite
    }

}

task alignPairedEnd {

    input {

        File input_fastq1
        File input_fastq2
        File input_reference
        String output_bam
        String aligner
        String path_aligner
        String read_fastq_command
        Int min_cov
        Int num_threads
        Boolean sort
        String tmp_folder
        Boolean overwrite

    }

    command {

        set -e pipefail

        snaptools align-paired-end \
            --input-reference=${input_reference} \
            --input-fastq1=${input_fastq1} \
            --input-fastq2=${input_fastq2} \
            --output-bam=${output_bam} \
            --aligner=${aligner} \
            --path-to-aligner=${path_aligner} \
            --read-fastq-command=${read_fastq_command} \
            --min-cov=${min_cov} \
            --num-threads=${num_threads} \
            --tmp-folder=${tmp_folder} \
            --overwrite=${overwrite} \
            --if-sort=${sort}

    }

    output {

        String path_output = "${output_bam}"
    }

    runtime {
        docker: "hisplan/snaptools:latest"
    }

}
