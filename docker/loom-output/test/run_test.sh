#!/bin/bash

docker build ../ -t tmp-test-image
docker run -it --mount type=bind,source="${PWD}/data/",target=/data\
       --mount type=bind,source="${PWD}/testout",target=/testout\
       --rm tmp-test-image:latest python3 create_loom_optimus.py\
            --empty_drops_file /data/empty_drops_result.csv\
            --cell_metrics /data/merged-cell-metrics.csv.gz\
            --gene_metrics /data/merged-gene-metrics.csv.gz\
            --cell_id /data/sparse_counts_row_index.npy\
            --gene_id /data/sparse_counts_col_index.npy\
            --count_matrix /data/sparse_counts.npz\
            --annotation_file /data/gencode.v27.primary_assembly.annotation.gtf.gz \
            --output_path_for_loom /testout/output\
            --sample_id "testsample"\
#            --verbose
