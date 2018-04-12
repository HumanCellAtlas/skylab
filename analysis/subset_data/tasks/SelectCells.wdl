task SelectCells {
  File input_bam
  Int num_small_cell_barcodes
  Int num_large_cell_barcodes
  Int? min_reads_for_large_cell_barcode
  Int? min_reads_for_small_cell_barcode
  String tag = "CB"

  Int? disk

  command <<<
    set -e

    java -Xms2g -jar /root/picard.jar GroupReadsByTag \
    I=${input_bam} \
    O=histogram \
    TAG=${tag}

    # write largest num_large_cells into file
    cat histogram | grep -v "#" | awk '{if($2>${default="50" min_reads_for_large_cell_barcode})printf("%s\t%s\n", $1, $2)}' | sort -k2 -n | tail -n+1 | tail -${num_large_cell_barcodes} | cut -f1 >> tags_to_filter

    # write smallest num_small_cells into file
    cat histogram | grep -v "#" | awk '{if($2>${default="100" min_reads_for_small_cell_barcode})printf("%s\t%s\n", $1, $2)}' | sort -k2 -n | tail -n+1 | head -${num_small_cell_barcodes} | cut -f1 >> tags_to_filter

    if [ `wc -l tags_to_filter | awk '{print $1}'` -lt "${num_small_cell_barcodes + num_large_cell_barcodes}" ]; then
      echo "Did not find the right number of total cell barcodes, expected ${num_small_cell_barcodes + num_large_cell_barcodes} but got `wc -l histogram | awk '{print $1}'`"
      exit 1
    fi

    if [ `uniq -d tags_to_filter | wc -l | awk '{print $1}'` -gt "1" ]; then
      echo "Not all cell barcodes in final file are unique, make sure there are enough cell barcodes with given parameters to select only unique barcodes"
      exit 1
    fi
  >>>
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    disks: "local-disk " + if defined(disk) then disk + " HDD" else ceil(size(input_bam, "GB")) + 10 + " HDD"
    memory: "3 GB"
  }
  output {
    Array[String] tag_values = read_lines("tags_to_filter")
    File histogram = "histogram"
  }
}
