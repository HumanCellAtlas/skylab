task FastQC {
    Array[File] fastq_files
    File? limits_file
    Int startRead = 250000
    Int nRead = 250000
    String docker = "quay.io/biocontainers/fastqc:0.11.8--1"
    Int machine_mem_mb = 3850
    Int disk = 100
    Int preemptible = 3

    String dollar = "$"
    parameter_meta {
        fastq_files: "input fastq files"
        limits_file: "(optional) limits file to use with fastqc"
        startRead: "(optional) start fastqc at the nth read of the file"
        nRead: "(optional) use (at most) n reads for fastqc"
        docker: "(optional) the docker image containing the runtime environment for this task"
        disk: "(optional) the amount of disk space (GiB) to provision for this task"
        preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
        machine_mem_mb: "(optional) the amount of memory (MiB) to provision for this task"

    }

    command <<<
        set -e

        mkdir outputs
        declare -a fastqs=()
        for fastq in ${sep=' ' fastq_files}
        do
            outname=`basename ${dollar}fastq .fastq.gz`_skip${startRead}_read${nRead}.fastq
            zcat ${dollar}fastq | head -n ${4*(startRead + nRead)} | tail -n ${4*nRead} > ${dollar}outname
            fastqs+=(${dollar}outname)
        done

        fastqc ${dollar}{fastqs[@]} -o outputs ${"--limits " + limits_file}
    >>>

    runtime {
        docker: docker
        memory: "${machine_mem_mb} MiB"
        disks: "local-disk ${disk} HDD"
        preemptible: preemptible
    }

    output {
        Array[File] fastqc_htmls = glob("outputs/*.html")
        Array[File] fastqc_zips = glob("outputs/*.zip")
    }
}