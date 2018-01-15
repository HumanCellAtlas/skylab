import "SelectCells.wdl" as selectcells
import "SubsetBamByTag.wdl" as subsetbambytag
import "SplitBam.wdl" as splitbam
import "SamToFastq.wdl" as samtofastq

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
  call selectcells.SelectCells as Select {
    input:
      input_bam = input_bam,
      num_large_cell_barcodes = num_large_cell_barcodes,
      num_small_cell_barcodes = num_small_cell_barcodes,
      min_reads_for_large_cell_barcode = min_reads_for_large_cell_barcode,
      min_reads_for_small_cell_barcode = min_reads_for_small_cell_barcode
  }

  # Grab reads from selected cell barcodes
  call subsetbambytag.SubsetBamByTag as SubsetBamTag {
    input:
      input_bam = input_bam,
      tag_values = Select.tag_values,
  }

  Int subset_bam_size = ceil(size(SubsetBamTag.subset_bam, "GB"))
  # Split bam into bams described by barcode_tags_to_split_by, quality_tags_to_split_by, separators_to_split_by
  call splitbam.SplitBam as Split {
    input:
      subset_bam = SubsetBamTag.subset_bam,
      barcode_tags_to_split_by = barcode_tags_to_split_by,
      quality_tags_to_split_by = quality_tags_to_split_by,
      separators_to_split_by = separators_to_split_by
  }

  # Revert bams into fastqs
  scatter (index in range(length(Split.split_bams))) {
    String base = if (defined(final_fastq_basenames)) then select_first([final_fastq_basenames])[index] else basename(Split.split_bams[index])
    Int split_bam_size = ceil(size(Split.split_bams[index], "GB"))

    call samtofastq.SamToFastq {
      input:
        input_bam = Split.split_bams[index],
        output_basename = base
    }
  }
}
