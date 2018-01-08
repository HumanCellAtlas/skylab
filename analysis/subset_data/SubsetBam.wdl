workflow SubsetBam {
  # Bam to subset
  File aligned_gold_bam

  # How many cells to grab from the largest and smallest end
  Int num_large_cells
  Int num_small_cells

  # How many reads a cell must have before it can be selected
  Int? min_reads_per_cell

  # Array of comma separated tag values to put in each fastq e.g ["ORIG", "SR", "CB,UB"] "ORIG" means to
  # grab the original sequence/qualities from the read
  Array[String] barcode_tags_to_split_by
  Array[String] quality_tags_to_split_by
  Array[String]? separators_to_split_by

  Array[String]? final_fastq_basenames

  Int bam_size = ceil(size(aligned_gold_bam, "GB"))
  Int additional_disk = 10

  # Figure out which cell barcodes to select from bam
  call SelectCells {
    input:
      aligned_bam = aligned_gold_bam,
      num_large_cells = num_large_cells,
      num_small_cells = num_small_cells,
      min_reads_per_cell = min_reads_per_cell,
      disk = bam_size + additional_disk
  }

  # Grab reads from selected cell barcodes
  call SubsetBamByTag {
    input:
      aligned_bam = aligned_gold_bam,
      tag_values = SelectCells.tag_values,
      disk = (bam_size * 2) + additional_disk # 2 * full bam size here as that would be the max subset size...
  }

  Int subset_bam_size = ceil(size(SubsetBamByTag.subset_bam, "GB"))
  # Split bam into bams described by barcode_tags_to_split_by, quality_tags_to_split_by, separators_to_split_by
  call SplitBam {
    input:
    subset_bam = SubsetBamByTag.subset_bam,
    barcode_tags_to_split_by = barcode_tags_to_split_by,
    quality_tags_to_split_by = quality_tags_to_split_by,
    separators_to_split_by = separators_to_split_by,
    disk = (subset_bam_size * 2) + additional_disk # max split size would be 2 * subset size
  }

  # Revert bams into fastqs
  scatter (index in range(length(SplitBam.split_bams))) {
    String base = if (defined(final_fastq_basenames)) then select_first([final_fastq_basenames])[index] else basename(SplitBam.split_bams[index])
    Int split_bam_size = ceil(size(SplitBam.split_bams[index], "GB"))

    call SamToFastq {
      input:
        bam = SplitBam.split_bams[index],
        output_basename = base,
        disk = (split_bam_size * 2) + additional_disk  # fast q shouldn't be larger than the bam size
    }
  }
}

task SelectCells {
  File aligned_bam
  Int num_large_cells
  Int num_small_cells
  Int? min_reads_per_cell
  String tag = "CB"

  Int disk

  command <<<
    set -e

    java -jar /root/picard.jar GroupReadsByTag \
    I=${aligned_bam} \
    O=histogram \
    T=${tag}

    # write largest num_large_cells into file
    cat histogram | grep -v "#" | awk '{if($2>${default="100" min_reads_per_cell})printf("%s\t%s\n", $1, $2)}' | sort -k2 -n | tail -n+1 | tail -${num_large_cells} | cut -f1 >> tags_to_filter

    # write smallest num_small_cells into file
    cat histogram | grep -v "#" | awk '{if($2>${default="100" min_reads_per_cell})printf("%s\t%s\n", $1, $2)}' | sort -k2 -n | tail -n+1 | head -${num_small_cells} | cut -f1 >> tags_to_filter
  >>>
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    disks: "local-disk " + disk + " HDD"
  }
  output {
    Array[String] tag_values = read_lines("tags_to_filter")
    File histogram = "histogram"
  }
}

task SubsetBamByTag {
  File aligned_bam
  Array[String] tag_values
  String tag = "CB"

  Int disk

  command {
    java -jar /root/picard.jar FilterSamReads \
    I=${aligned_bam} \
    FILTER=includeTagValues \
    O=selected_cells.bam \
    T=${tag} \
    TV=${sep=" TV=" tag_values}
  }
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    disks: "local-disk " + disk + " HDD"
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

  Int disk

  command {
    java -Xmx3g -jar /root/picard.jar SplitBamByTags \
    T=${sep=" T=" barcode_tags_to_split_by} \
    Q=${sep=" Q=" quality_tags_to_split_by} \
    I=${subset_bam} \
    O=split \
    ${separator}${default="" sep=" S=" separators_to_split_by}
  }
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    disks: "local-disk " + disk + " HDD"
    memory: "7 GB"
  }
  output {
    Array[File] split_bams = glob("split*")
  }
}

task SamToFastq {
  File bam
  String output_basename

  Int disk

  command {
    java -jar /root/picard.jar SamToFastq \
    I=${bam} \
    F=${output_basename}.fastq
  }
  runtime {
    docker: "jsotobroad/gold_data:1.0"
    disks: "local-disk " + disk + " HDD"
  }
  output {
    File fastq = "${output_basename}.fastq"
  }
}
