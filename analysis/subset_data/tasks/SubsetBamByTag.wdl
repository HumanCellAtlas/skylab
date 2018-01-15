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
