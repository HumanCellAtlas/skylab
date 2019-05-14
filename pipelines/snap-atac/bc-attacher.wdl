version 1.0

workflow BCAttacher {
    input {
        File index1
        File index2
        File read1
        File read2
    }

    parameter_meta {
        index1: "index 1 input fastq"
        index2: "index 2 input fastq"
        read1: "read 1 input fastq"
        read2: "read 2 input fastq"
    }

    call AttachBarcode {
        input:
            index1 = index1,
            index2 = index2,
            read1 = read1,
            read2 = read2
    }

    output {
        File output_read1 = AttachBarcode.output_read1
        File output_read2 = AttachBarcode.output_read2
    }
}

task AttachBarcode {
    input {
        File index1
        File index2
        File read1
        File read2
        String out_read1_filename = "R1.fastq.gz"
        String out_read2_filename = "R2.fastq.gz"
        String docker_image = "hisplan/bc-attacher"
    }

    Int num_threads = 1
    Float input_size = size(index1, "GiB") + size(index2, "GiB") + size(read1, "GiB") + size(read2, "GiB")

    command {
        set -euo pipefail

        python bc_attacher.py \
            --index1 ~{index1} \
            --index2 ~{index2} \
            --read1 ~{read1} \
            --read2 ~{read2} \
            --out-read1 ~{out_read1_filename} \
            --out-read2 ~{out_read2_filename}
    }

    output {
        File output_read1 = out_read1_filename
        File output_read2 = out_read2_filename
    }

    runtime {
        docker: docker_image
        cpu: num_threads
        memory: "16 GB"
        disks: "local-disk " + ceil(3.5 * (if input_size < 1 then 1 else input_size )) + " HDD"
    }
}
