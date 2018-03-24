import "FastqToUBam.wdl" as fq2bam
import "Attach10xBarcodes.wdl" as attach
import "SplitBamByCellBarcode.wdl" as split
import "MergeSortBam.wdl" as merge
import "CreateCountMatrix.wdl" as count
import "StarAlignBamSingleEnd.wdl" as StarAlignBamWDL
import "TagGeneExon.wdl" as TagGeneExonWDL
import "CorrectUmiMarkDuplicates.wdl" as CorrectUmiMarkDuplicatesWDL
import "SequenceDataWithMoleculeTagMetrics.wdl" as metrics
import "TagSortBam.wdl" as TagSortBam

# The optimus 3' pipeline processes 10x genomics sequencing data based on the v2
# chemistry. It corrects cell barcodes and UMIs, aligns reads, marks duplicates, and
# returns data as alignments in BAM format and as counts in sparse matrix exchange
# format.
workflow Optimus {

  # Sequencing data inputs (fastq)
  Array[File] r1 # forward read, contains cell barcodes and molecule barcodes
  Array[File] r2 # reverse read, contains cDNA fragment generated from captured mRNA
  Array[File]? i1 # (optional) index read, for demultiplexing of multiple samples on one flow cell.

  String sample_id  # name of sample matching this file, inserted into read group header

  # organism reference parameters
  File tar_star_reference  # star reference
  File annotations_gtf  # gtf containing annotations for gene tagging
  File ref_genome_fasta  # genome fasta file

  # 10x v2 parameters
  File whitelist  # 10x genomics cell barcode whitelist for 10x V2

  # environment-specific parameters
  String fastq_suffix = ""  # when running in green box, need to add ".gz" for picard to detect
  Array[Int] indices = range(length(r1)) # this scatters matched [r1, r2, i1] fastq arrays

  scatter (index in indices) {
    call fq2bam.FastqToUBam {
      input:
        fastq_file = r2[index],
        sample_id = sample_id,
        fastq_suffix = fastq_suffix
    }

    # if the index is passed, attach it to the bam file
    if (defined(i1)) {
      Array[File] non_optional_i1 = select_first([i1])
      call attach.Attach10xBarcodes as attach_barcodes {
        input:
          r1 = r1[index],
          i1 = non_optional_i1[index],
          u2 = FastqToUBam.bam_output,
          whitelist = whitelist
      }
    }

    # if the index is not passed, proceed without it.
    if (!defined(i1)) {
      call attach.Attach10xBarcodes as attach_barcodes_no_index {
        input:
          r1 = r1[index],
          u2 = FastqToUBam.bam_output,
          whitelist = whitelist
      }
    }

    File barcoded_bam = select_first([attach_barcodes.bam_output, attach_barcodes_no_index.bam_output])
  }

  call merge.MergeSortBamFiles as MergeUnsorted {
    input:
      bam_inputs = barcoded_bam,
      sort_order = "unsorted"
  }

  call split.SplitBamByCellBarcode {
    input:
      bam_input = MergeUnsorted.output_bam
  }

  scatter (bam in SplitBamByCellBarcode.bam_output_array) {
    call StarAlignBamWDL.StarAlignBamSingleEnd as StarAlign {
      input:
        bam_input = bam,
        tar_star_reference = tar_star_reference
    }

    call TagGeneExonWDL.TagGeneExon as TagGenes {
      input:
        bam_input = StarAlign.bam_output,
        annotations_gtf = annotations_gtf
    }

    call CorrectUmiMarkDuplicatesWDL.CorrectUmiMarkDuplicates as CorrectUmiMarkDuplicates {
      input:
        bam_input = TagGenes.bam_output
    }

    call TagSortBam.GeneSortBam {
      input:
        bam_input = CorrectUmiMarkDuplicates.bam_output
    }

    call TagSortBam.CellSortBam {
      input:
        bam_input = CorrectUmiMarkDuplicates.bam_output
    }

    call metrics.CalculateGeneMetrics {
      input:
        bam_input = GeneSortBam.bam_output
    }

    call metrics.CalculateCellMetrics {
      input:
        bam_input = CellSortBam.bam_output
    }
  }

  call merge.MergeSortBamFiles as MergeSorted {
    input:
      bam_inputs = CorrectUmiMarkDuplicates.bam_output,
      sort_order = "coordinate"
  }

  call count.DropSeqToolsDigitalExpression {
    input:
      bam_input = MergeSorted.output_bam,
      whitelist = whitelist
  }

  call metrics.MergeGeneMetrics {
    input:
      metric_files = CalculateGeneMetrics.gene_metrics
  }

  call metrics.MergeCellMetrics {
    input:
      metric_files = CalculateCellMetrics.cell_metrics
  }

  output {
      File bam = MergeSorted.output_bam
      File matrix = DropSeqToolsDigitalExpression.matrix_output
      File matrix_summary = DropSeqToolsDigitalExpression.matrix_summary
      File cell_metrics = MergeCellMetrics.cell_metrics
      File gene_metrics = MergeGeneMetrics.gene_metrics
  }
}
