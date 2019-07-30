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

     call ValidateBam as ValidateBam {
         input:
            bam = bam,
            expected_checksum = expected_bam_hash
     }

     call ValidateMatrix as ValidateMatrix {
         input:
             matrix = matrix,
             matrix_row_index = matrix_row_index,
	     matrix_col_index = matrix_col_index,
             expected_matrix_hash = expected_matrix_hash
     }

     call ValidateMetrics {
         input:
             cell_metrics = cell_metrics,
             gene_metrics = gene_metrics,
             
             expected_cell_metric_hash = expected_cell_metric_hash,
             expected_gene_metric_hash = expected_gene_metric_hash
     }

     call GenerateReport as GenerateReport {
         input:
             bam_validation_result = ValidateBam.result,
             matrix_validation_result = ValidateMatrix.result,
             metric_and_index_validation_result = ValidateMetrics.result
    }
}

task ValidateBam {
    File bam
    String expected_checksum
    Int required_disk = ceil(size(bam, "G") * 1.1)

    command <<<
        echo Starting checksum generation...
        # calculate hash for alignment positions only (a reduced bam hash)
        samtools view -F 256 "${bam}" | cut -f 1-11 | md5sum | awk '{print $1}' > md5_checksum_reduced.txt
        echo Reduced checksum generation complete

        calculated_checksum = $(cat md5_checksum_reduced.txt)

        if [ "$calculated_checksum" == ${expected_checksum} ]
        then
             echo Computed and expected bam hashes match ( $calculated_checksum )
             printf PASS > result.txt
        else 
             echo Computed ( $calculated_checksum ) and expected ( ${expected_checksum} ) hashes do not match
             printf FAIL > result.txt
        fi
    >>>
  
    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string("result.txt")
    }
}

task ValidateMatrix {
    File matrix
    File matrix_row_index
    File matrix_col_index
    String expected_matrix_hash
    
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
       computed_hash=$(find . -name "*.npy" -type f -exec md5sum {} \; | sort -k 2 | md5sum | awk '{print $1}')
       cd ..

       if [ $computed_hash == ${expected_matrix_hash} ]; then
           echo Computed and expected matrix hashes match ( $computed_hash )
           printf PASS > result.txt
       else 
           echo Computed hash ( $computed_hash ) did not match expected matrix hash ( ${expected_matrix_hash} )
           printf FAIL > result.txt
       fi

    >>>
  
    runtime {
        docker: "quay.io/humancellatlas/secondary-analysis-samtools:v0.2.2-1.6"
        cpu: 1
        memory: "3.75 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string('result.txt')
    }  

}

task ValidateMetrics {
    File cell_metrics
    File gene_metrics

    String expected_cell_metric_hash
    String expected_gene_metric_hash

    Int required_disk = ceil( (size(cell_metrics, "G") + size(gene_metrics, "G") )* 1.1)

    command <<<   
       set -eo pipefail

        # check matrix row and column indexes files hash
        gene_metric_hash=$(zcat "${gene_metrics}" | md5sum | awk '{print $1}' > gene_metric_hash.txt)
        cell_metric_hash=$(zcat "${cell_metrics}" | md5sum | awk '{print $1}' > cell_metric_hash.txt)

        fail=false

        if [ $gene_metric_hash != ${expected_gene_metric_hash} ]; then
            echo Computed ( $gene_metric_hash ) and expected ( ${expected_gene_metric_hash} ) checksums do not match
            fail=true
        else 
            echo Computed and expected gene metrics match ( $gene_metric_hash )
        fi

        if [ $cell_metric_hash != ${expected_cell_metric_hash} ]; then
            echo Computed ( $cell_metric_hash ) and expected ( $expected_cell_metric_hash ) cell metrics hashes don't match
            fail=true
        else 
            echo Computed and expected cell metrics match ( $cell_metric_hash )
        fi

        if [ $fail ]; then
            printf FAIL > result.txt
        else 
            printf PASS > result.txt
        fi

    >>>

    runtime {
        docker: "ubuntu:18.04"
        cpu: 1
        memory: "1.00 GB"
        disks: "local-disk ${required_disk} HDD"
    }

    output {
        String result = read_string("result.txt")
    }

}

task GenerateReport {
      String bam_validation_result
      String metric_and_index_validation_result
      String matrix_validation_result

      Int required_disk = 1

  command <<<

    set -eo pipefail

    # test each output for equality, echoing any failure states to stdout
    fail=false

    echo Bam Validation: ${bam_validation_result}
    if [ ${bam_validation_result} == "FAIL" ]; then
        fail=true
    fi

    echo Metrics and indices Validation: ${metric_and_index_validation_result}
    if [ ${metric_and_index_validation_result} == "FAIL" ]; then
        fail=true
    fi
    
    echo Matrix Validation: ${matrix_validation_result}
    if [ ${matrix_validation_result} == "FAIL" ]; then
        fail=true
    fi

    if [ $fail == "true" ]; then exit 1; fi

  >>>
  
  runtime {
    docker: "ubuntu:18.04"
    cpu: 1
    memory: "1.0 GB"
    disks: "local-disk ${required_disk} HDD"
  }
}
