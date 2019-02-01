task MarkDuplicatesUmiTools {
    File bam_input

    # runtime values
    String docker = "quay.io/humancellatlas/secondary-analysis-umitools:0.0.1"

    ## TODO: Optimize these values
    Int machine_mem_mb = 10000
    Int cpu = 1
    Int disk = ceil(size(bam_input, "G") * 6) + 50
    Int preemptible = 0

    meta {
        description: "Marks duplicates using umitools group specifically for single-cell experiments"
    }

    command {
        set -e

        umi_tools group \
            -I ${bam_input} \
            -L outlog.txt \
            -E outerr.txt \
            -S duplicate_marked.bam \
            --output-bam \
            --extract-umi-method=tag \
            --umi-tag UR \
            --method directional \
            --per-gene \
            --per-cell \
            --cell-tag CB \
            --gene-tag GE
            --no-sort-output \
            --group-out groupout.txt \
            --umi-group-tag UB



}