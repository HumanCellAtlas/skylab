task ScatterBam {

    File bam_to_scatter
    Int scatter_width

    Int disk_size = ceil(size(bam_to_scatter, "GB") * 3)

    command <<<
        mkdir splitted_bams
        java -Xms7g -jar /usr/gitc/picard.jar \
            SplitSamByNumberOfReads \
            INPUT=${bam_to_scatter} \
            SPLIT_TO_N_FILES=${scatter_width} \
            OUTPUT_PREFIX=$(basename ${bam_to_scatter} .bam).split \
            OUTPUT=splitted_bams
    >>>

    output {
        Array[File] splitted_bams = glob("splitted_bams/*")
    }

    runtime {
        disk: "local-disk ${disk_size} HDD"
        cpu: 2
        memory: "7.5 GB"
        docker: "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.3.3-1513176735"
    }
}