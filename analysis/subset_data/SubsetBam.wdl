workflow SubsetBam {
  # Bam to subset
  File input_bam

  # How many cells to grab from the largest and smallest end
  Int num_large_cell_barcodes
  Int num_small_cell_barcodes

  # How many reads a cell must have before it can be selected
  Int? min_reads_for_small_cell_barcode
  Int? min_reads_for_large_cell_barcode

  # Array of comma separated tag values to put in each fastq e.g ["ORIG", "SR", "CB,UB"] "ORIG" means to
  # grab the original sequence/qualities from the read
  # the length of these arrays should all be equal
  Array[String] barcode_tags_to_split_by
  Array[String] quality_tags_to_split_by
  Array[String]? separators_to_split_by

  Array[String]? final_fastq_basenames

  # Figure out which cell barcodes to select from bam
  call SelectCells {
    input:
      input_bam = input_bam,
      num_large_cell_barcodes = num_large_cell_barcodes,
      num_small_cell_barcodes = num_small_cell_barcodes,
      min_reads_for_large_cell_barcode = min_reads_for_large_cell_barcode,
      min_reads_for_small_cell_barcode = min_reads_for_small_cell_barcode
  }

  # Grab reads from selected cell barcodes
  call SubsetBamByTag {
    input:
      input_bam = input_bam,
      tag_values = SelectCells.tag_values,
  }

  Int subset_bam_size = ceil(size(SubsetBamByTag.subset_bam, "GB"))
  # Split bam into bams described by barcode_tags_to_split_by, quality_tags_to_split_by, separators_to_split_by
  call SplitBam {
    input:
      subset_bam = SubsetBamByTag.subset_bam,
      barcode_tags_to_split_by = barcode_tags_to_split_by,
      quality_tags_to_split_by = quality_tags_to_split_by,
      separators_to_split_by = separators_to_split_by
  }

  # Revert bams into fastqs
  scatter (index in range(length(SplitBam.split_bams))) {
    String base = if (defined(final_fastq_basenames)) then select_first([final_fastq_basenames])[index] else basename(SplitBam.split_bams[index])
    Int split_bam_size = ceil(size(SplitBam.split_bams[index], "GB"))

    call SamToFastq {
      input:
        input_bam = SplitBam.split_bams[index],
        output_basename = base
    }
  }
}

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

task SubsetBamByTag {
  File input_bam
  Array[String] tag_values
  String tag = "CB"

  Int? disk

  command {
    java -Xms2g -jar /root/picard.jar FilterSamReads \
    I=${input_bam} \
    FILTER=includeTagValues \
    O=selected_cells.bam \
    TAG=${tag} \
    TAG_VALUE=${sep=" TAG_VALUE=" tag_values}
  }
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    # 2 * full bam size here as that would be the max subset size
    disks: "local-disk " + if defined(disk) then disk + " HDD" else (ceil(size(input_bam, "GB")) * 2) + 10 + " HDD"
    memory: "3 GB"
  }
  output {
    File subset_bam = "selected_cells.bam"
  }
}


task SplitBam {
  File subset_bam
  Array[String] barcode_tags_to_split_by
  Array[String] quality_tags_to_split_by
  Array[String]? separators_to_split_by
  String separator = if defined(separators_to_split_by) then "S= " else ""

  Int? disk

  command {
    # T - comma separated list of tags with base sequence values used to create bams with those sequences.
    # the value of "ORIG" will use the original sequence from the read
    #  e.g. ["CR,UR", "ORIG"] will create two bams, one whos sequence is a concatenation of the values of the CR and UR tags
    # Q - equivalent to T but instead of base sequence values this is for base quality values.  The length of this array
    # must be the same length as T
    # S - optional base sequence value to separate every set of tags pass for each split bam.  The length of this array
    # must be the same length as T
    java -Xms6g -jar /root/picard.jar SplitBamByTags \
    T=${sep=" T=" barcode_tags_to_split_by} \
    Q=${sep=" Q=" quality_tags_to_split_by} \
    I=${subset_bam} \
    O=split \
    ${separator}${default="" sep=" S=" separators_to_split_by}
  }
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    # max split size would be 2 * subset size
    disks: "local-disk " + if defined(disk) then disk + " HDD" else (ceil(size(subset_bam, "GB")) * 2) + 10 + " HDD"
    memory: "7 GB"
  }
  output {
    Array[File] split_bams = glob("split*")
  }
}

task SamToFastq {
  File input_bam
  String output_basename

  Int? disk

  command {
    java -Xms2g -jar /root/picard.jar SamToFastq \
    I=${input_bam} \
    F=${output_basename}.fastq
  }
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    # fast q shouldn't be larger than the bam size
    disks: "local-disk " + if defined(disk) then disk + " HDD" else (ceil(size(input_bam, "GB")) * 2) + 10 + " HDD"
    memory: "3 GB"
  }
  output {
    File fastq = "${output_basename}.fastq"
  }
}
