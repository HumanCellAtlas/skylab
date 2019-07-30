workflow ValidateOptimus {
     meta {
         description: "Validate Optimus Outputs"
     }

     File bam
     File matrix
     File matrix_row_index
     File matrix_col_index
     File cell_metrics
     File gene_metrics

     String expected_bam_hash
     String expected_matrix_hash
     String expected_matrix_row_hash
     String expected_matrix_col_hash
     String expected_cell_metric_hash
     String expected_gene_metric_hash
     String expected_reduced_bam_hash

     call ValidateBam as ValidateBam {
         input:
            bam = bam
     }

     call ValidateMatrix as ValidateMatrix {
         input:
             matrix = matrix
     }

     call ValidateMetricsAndIndexes {
         input:
             matrix_row_index = matrix_row_index,
             matrix_col_index = matrix_col_index,
             cell_metrics = cell_metrics,
             gene_metrics = gene_metrics
     }

     call MainOptimusValidation as MainOptimusValidation {
         input:
             bam_hash = ValidateBam.md5sum,
             bam_reduced_hash = ValidateBam.md5sum_reduced,
             matrix_hash = ValidateMatrix.matrix_hash,

	     matrix_row_index_hash = ValidateMetricsAndIndexes.matrix_row_index_hash,
             matrix_col_index_hash = ValidateMetricsAndIndexes.matrix_col_index_hash,
             cell_metrics_hash = ValidateMetricsAndIndexes.cell_metrics_hash,
             gene_metrics_hash = ValidateMetricsAndIndexes.gene_metrics_hash,

             expected_bam_hash = expected_bam_hash,
             expected_matrix_hash = expected_matrix_hash,
             expected_matrix_row_hash = expected_matrix_row_hash,
             expected_matrix_col_hash = expected_matrix_col_hash,
             expected_cell_metric_hash = expected_cell_metric_hash,
             expected_gene_metric_hash = expected_gene_metric_hash,
             expected_reduced_bam_hash = expected_reduced_bam_hash
    }
}

task ValidateBam {
    File bam
    Int required_disk = ceil(size(bam, "G") * 1.1)

    command <<<
        # calculate hash as above, but ignore run-specific bam headers
        samtools view "${bam}" | md5sum | awk '{print $1}' > md5_checksum.txt
    
        # calculate hash for alignment positions only (a reduced bam hash)
        samtools view -F 256 "${bam}" | cut -f 1-11 | md5sum | awk '{print $1}' > md5_checksum_reduced.txt
    >>>
  
    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String md5sum = read_string("md5_checksum.txt")
        String md5sum_reduced = read_string("md5_checksum_reduced.txt")
    }
}

task ValidateMatrix {
    File matrix
    
    Int required_disk = ceil( size(matrix, "G") * 1.1 )

    command <<<
       set -eo pipefail

       ## Generate md5sum for the matrix.npz file, the file is first deflated
       ## to avoid changes in the compression process from affecting output
       ## the process is performed in new directory to avoid other *.npy files
       ## that may be inputs from affecting output
       mkdir matrix_deflate
       cp "${matrix}" matrix_deflate/matrix.npz
       cd matrix_deflate
       unzip matrix.npz
       rm matrix.npz
       find . -name "*.npy" -type f -exec md5sum {} \; | sort -k 2 | md5sum | awk '{print $1}' > ../matrix_hash.txt
       cd ..

    >>>
  
    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String matrix_hash = read_string("matrix_hash.txt")
    }  

}

task ValidateMetricsAndIndexes {
    File matrix_row_index
    File matrix_col_index
    File cell_metrics
    File gene_metrics

    Int required_disk = ceil((size(matrix_row_index, "G") + size(matrix_col_index, "G") + size(cell_metrics, "G") + size(gene_metrics, "G") )* 1.1)

    command <<<   
       set -eo pipefail

        # check matrix row and column indexes files hash
        cat "${matrix_row_index}" | md5sum | awk '{print $1} > matrix_row_index_hash.txt
        cat "${matrix_col_index}" | md5sum | awk '{print $1} > matrix_col_index_hash.txt
    
        # check gene and cell metrics
        zcat "${gene_metrics}" | md5sum | awk '{print $1} > gene_metric_hash.txt
        zcat "${cell_metrics}" | md5sum | awk '{print $1} > cell_metric_hash.txt
    >>>

    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String matrix_row_index_hash = read_string("matrix_row_index.txt")
        String matrix_col_index_hash = read_string("matrix_col_index.txt")
        String gene_metrics_hash = read_string("gene_metric_hash.txt")
        String cell_metrics_hash = read_string("cell_metric_hash.txt")
    }

}

task MainOptimusValidation {
      String bam_hash
      String bam_reduced_hash
      String matrix_hash

      String matrix_row_index_hash
      String matrix_col_index_hash
      String gene_metrics_hash
      String cell_metrics_hash

      Int required_disk = 2

      String expected_bam_hash
      String expected_matrix_hash
      String expected_matrix_row_hash
      String expected_matrix_col_hash
      String expected_cell_metric_hash
      String expected_gene_metric_hash
      String expected_reduced_bam_hash

  command <<<

    set -eo pipefail

    # test each output for equality, echoing any failure states to stdout
    fail=false

    if [ "$bam_reduced_hash" != "${expected_reduced_bam_hash}" ]; then
      >&2 echo "bam_reduced_hash ($bam_reduced_hash) did not match expected hash (${expected_reduced_bam_hash})"
      fail=true
    fi

    if [ "$matrix_hash" != "${expected_matrix_hash}" ]; then
      >&2 echo "matrix_hash ($matrix_hash) did not match expected hash (${expected_matrix_hash})"
      fail=true
    fi

    if [ "$matrix_row_index_hash" != "${expected_matrix_row_hash}" ]; then
      >&2 echo "matrix_row_index_hash ($matrix_row_index_hash) did not match expected hash (${expected_matrix_row_hash})"
      fail=true
    fi

    if [ "$matrix_col_index_hash" != "${expected_matrix_col_hash}" ]; then
      >&2 echo "matrix_col_index_hash ($matrix_col_index_hash) did not match expected hash (${expected_matrix_col_hash})"
      fail=true
    fi

    if [ "$bam_hash" != "${expected_bam_hash}" ]; then
      >&2 echo "bam_hash ($bam_hash) did not match expected hash (${expected_bam_hash})"
      fail=true
    fi

    if [ "$gene_metric_hash" != "${expected_gene_metric_hash}" ]; then
      >&2 echo "gene_metric_hash ($gene_metric_hash) did not match expected hash (${expected_gene_metric_hash})"
      fail=true
    fi

    if [ "$cell_metric_hash" != "${expected_cell_metric_hash}" ]; then
      >&2 echo "cell_metric_hash ($cell_metric_hash) did not match expected hash (${expected_cell_metric_hash})"
      fail=true
    fi

    if [ $fail == "true" ]; then exit 1; fi

  >>>
  
  runtime {
    docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
    cpu: 1
    memory: "3.75 GB"
    disks: "local-disk ${required_disk} HDD"
  }
}
