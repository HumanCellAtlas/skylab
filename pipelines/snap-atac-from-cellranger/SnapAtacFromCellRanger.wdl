version 1.0

workflow SnapAtacFromCellRanger {
    input {
        File possorted_bam
        String genome_name
        File genome_size_file
    }

    parameter_meta {
        possorted_bam: "position-sorted BAM file produced by Cell Ranger ATAC"
        genome_name: "genome name (e.g. hg38, mm10)"
        genome_size_file: "text file that contains corresponding genome sizes (e.g. hg38.chrom.sizes)"
    }

    call GetPrefix {
        input:
            bam_filename = possorted_bam
    }

    call AttachBarcodes {
        input:
            possorted_bam = possorted_bam,
            prefix = GetPrefix.prefix
    }

    call SnapPreprocessing {
        input:
            snap_sorted_bam = AttachBarcodes.snap_sorted_bam,
            prefix = GetPrefix.prefix,
            genome_name = genome_name,
            genome_size_file = genome_size_file
    }

    call AddCellByBinMatrix {
        input:
            snap_file = SnapPreprocessing.snap_file
    }

    output {
        File snap_sorted_bam = AttachBarcodes.snap_sorted_bam
        File snap_file = AddCellByBinMatrix.out_snap_file
    }
}

task GetPrefix {
    input {
        String bam_filename
        String docker_image = "python:3.6"
    }

    command {
        set -euo pipefail

        python <<CODE
        import os
        print(os.path.splitext(os.path.basename("~{bam_filename}"))[0])
        CODE
    }

    output {
        String prefix = read_string(stdout())
    }

    runtime {
        docker: docker_image
    }
}

task AttachBarcodes {
    input {
        File possorted_bam
        String prefix
        String docker_image = "quay.io/broadinstitute/samtools:1.9"
    }

    Int num_threads = 5
    Int ram_size_gb = 1

    Float input_size = size(possorted_bam, "GB")

    String header_filename = prefix + ".header"
    String snap_bam_filename = prefix + ".snap.bam"
    String snap_sorted_bam_filename = prefix + ".snap.sorted.bam"

    command <<<
        set -euo pipefail

        # extract the header
        samtools view ~{possorted_bam} -H > ~{header_filename}

        # create a bam file with the barcode embedded into the read name
        # e.g.
        # GGTTGCGAGCCGCAAA-1:A00519:218:HJYFLDSXX:2:1216:26458:34976
        cat <( cat ~{header_filename} ) \
            <( samtools view ~{possorted_bam} | awk '{for (i=12; i<=NF; ++i) { if ($i ~ "^CB:Z:"){ td[substr($i,1,2)] = substr($i,6,length($i)-5); } }; printf "%s:%s\n", td["CB"], $0 }' ) \
            | samtools view -bS - > ~{snap_bam_filename}

        samtools sort -n -@ ~{num_threads} -m ~{ram_size_gb}G \
            ~{snap_bam_filename} \
            -o ~{snap_sorted_bam_filename}
    >>>

    output {
        File snap_sorted_bam = snap_sorted_bam_filename
    }

    runtime {
        docker: docker_image
        # fixme:
        # disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
        disks: "local-disk 500 HDD"
        cpu: num_threads
        memory: "3.5 GB"
    }
}

task SnapPreprocessing {
    input {
        File snap_sorted_bam
        String prefix
        String genome_name
        File genome_size_file
        String docker_image = "hisplan/snaptools"
    }

    Int num_threads = 5

    Float input_size = size(snap_sorted_bam, "GB")

    String snap_filename = prefix + ".snap"

    command {
        set -euo pipefail

        snaptools snap-pre \
            --input-file=~{snap_sorted_bam} \
            --output-snap=~{snap_filename} \
            --genome-name=~{genome_name} \
            --genome-size=~{genome_size_file} \
            --min-mapq=30 \
            --min-flen=50 \
            --max-flen=1000 \
            --keep-chrm=TRUE \
            --keep-single=FALSE \
            --keep-secondary=False \
            --overwrite=True \
            --max-num=20000 \
            --min-cov=500 \
            --verbose=True
    }

    output {
        File snap_file = snap_filename
    }

    runtime {
        docker: docker_image
        # fixme:
        # disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
        disks: "local-disk 500 HDD"
        cpu: num_threads
        memory: "8 GB"
    }
}

task AddCellByBinMatrix {
    input {
        File snap_file
        String docker_image = "hisplan/snaptools"
    }

    Int num_threads = 5

    Float input_size = size(snap_file, "GB")

    command {
        set -euo pipefail

        snaptools snap-add-bmat	\
            --snap-file=~{snap_file} \
            --bin-size-lis 5000 100000
    }

    output {
        File out_snap_file = snap_file
    }

    runtime {
        docker: docker_image
        # fixme:
        # disks: "local-disk " + ceil(2 * (if input_size < 1 then 1 else input_size)) + " HDD"
        disks: "local-disk 500 HDD"
        cpu: num_threads
        memory: "8 GB"
    }
}
