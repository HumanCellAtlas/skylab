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
