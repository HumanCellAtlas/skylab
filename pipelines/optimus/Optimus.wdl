import "FastqToUBam.wdl" as FastqToUBam
import "Attach10xBarcodes.wdl" as Attach
import "SplitBamByCellBarcode.wdl" as Split
import "MergeSortBam.wdl" as Merge
import "CreateCountMatrix.wdl" as Count
import "StarAlignBamSingleEnd.wdl" as StarAlignBam
import "TagGeneExon.wdl" as TagGeneExon
import "CorrectUmiMarkDuplicates.wdl" as CorrectUmiMarkDuplicates
import "SequenceDataWithMoleculeTagMetrics.wdl" as Metrics
import "TagSortBam.wdl" as TagSortBam

workflow Optimus {
  meta {
    description: "The optimus 3' pipeline processes 10x genomics sequencing data based on the v2 chemistry. It corrects cell barcodes and UMIs, aligns reads, marks duplicates, and returns data as alignments in BAM format and as counts in sparse matrix exchange format."
  }

  # Sequencing data inputs
  Array[File] r1_fastq
  Array[File] r2_fastq
  Array[File]? i1_fastq
  String sample_id

  # organism reference parameters
  File tar_star_reference
  File annotations_gtf
  File ref_genome_fasta

  # 10x v2 parameters
  File whitelist

  # environment-specific parameters
  String fastq_suffix = ""
  # this is used to scatter matched [r1_fastq, r2_fastq, i1_fastq] arrays
  Array[Int] indices = range(length(r1_fastq))

  parameter_meta {
    r1_fastq: "forward read, contains cell barcodes and molecule barcodes"
    r2_fastq: "reverse read, contains cDNA fragment generated from captured mRNA"
    i1_fastq: "(optional) index read, for demultiplexing of multiple samples on one flow cell."
    sample_id: "name of sample matching this file, inserted into read group header"
    tar_star_reference: "star genome reference"
    annotations_gtf: "gtf containing annotations for gene tagging (must match star reference)"
    ref_genome_fasta: "genome fasta file (must match star reference)"
    whitelist: "10x genomics cell barcode whitelist for 10x V2"
    fastq_suffix: "when running in green box, need to add '.gz' for picard to detect the compression"
  }

  scatter (index in indices) {
    call FastqToUBam.FastqToUBam {
      input:
        fastq_file = r2_fastq[index],
        sample_id = sample_id,
        fastq_suffix = fastq_suffix
    }

    # if the index is passed, attach it to the bam file
    if (defined(i1_fastq)) {
      Array[File] non_optional_i1_fastq = select_first([i1_fastq])
      call Attach.Attach10xBarcodes as AttachBarcodes {
        input:
          r1_fastq = r1_fastq[index],
          i1_fastq = non_optional_i1_fastq[index],
          r2_unmapped_bam = FastqToUBam.bam_output,
          whitelist = whitelist
      }
    }

    # if the index is not passed, proceed without it.
    if (!defined(i1_fastq)) {
      call Attach.Attach10xBarcodes as AttachBarcodesNoIndex {
        input:
          r1_fastq = r1_fastq[index],
          r2_unmapped_bam = FastqToUBam.bam_output,
          whitelist = whitelist
      }
    }

    File barcoded_bam = select_first([AttachBarcodes.bam_output, AttachBarcodesNoIndex.bam_output])
  }

  call Merge.MergeSortBamFiles as MergeUnsorted {
    input:
      bam_inputs = barcoded_bam,
      sort_order = "unsorted"
  }

  call Split.SplitBamByCellBarcode {
    input:
      bam_input = MergeUnsorted.output_bam
  }

  scatter (bam in SplitBamByCellBarcode.bam_output_array) {
    call StarAlignBam.StarAlignBamSingleEnd as StarAlign {
      input:
        bam_input = bam,
        tar_star_reference = tar_star_reference
    }

    call TagGeneExon.TagGeneExon as TagGenes {
      input:
        bam_input = StarAlign.bam_output,
        annotations_gtf = annotations_gtf
    }

    call CorrectUmiMarkDuplicates.SortAndCorrectUmiMarkDuplicates {
      input:
        bam_input = TagGenes.bam_output
    }

    call TagSortBam.GeneSortBam {
      input:
        bam_input = SortAndCorrectUmiMarkDuplicates.bam_output
    }

    call TagSortBam.CellSortBam {
      input:
        bam_input = SortAndCorrectUmiMarkDuplicates.bam_output
    }

    call Metrics.CalculateGeneMetrics {
      input:
        bam_input = GeneSortBam.bam_output
    }

    call Metrics.CalculateCellMetrics {
      input:
        bam_input = CellSortBam.bam_output
    }

    call Count.CreateSparseCountMatrix {
      input:
        bam_input = CellSortBam.bam_output,
        gtf_file = annotations_gtf
    }
  }

  call Merge.MergeSortBamFiles as MergeSorted {
    input:
      bam_inputs = SortAndCorrectUmiMarkDuplicates.bam_output,
      sort_order = "coordinate"
  }

  call Metrics.MergeGeneMetrics {
    input:
      metric_files = CalculateGeneMetrics.gene_metrics
  }

  call Metrics.MergeCellMetrics {
    input:
      metric_files = CalculateCellMetrics.cell_metrics
  }

  call Count.MergeCountFiles {
    input:
      sparse_count_matrices = CreateSparseCountMatrix.sparse_count_matrix,
      row_indices = CreateSparseCountMatrix.row_index,
      col_indices = CreateSparseCountMatrix.col_index
  }

  output {
      File bam = MergeSorted.output_bam
      File matrix = MergeCountFiles.sparse_count_matrix
      File matrix_row_index = MergeCountFiles.row_index
      File matrix_col_index = MergeCountFiles.col_index
      File cell_metrics = MergeCellMetrics.cell_metrics
      File gene_metrics = MergeGeneMetrics.gene_metrics
  }
}
