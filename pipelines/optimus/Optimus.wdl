import "FastqToUBam.wdl" as FastqToUBam
import "Attach10xBarcodes.wdl" as Attach
import "SplitBamByCellBarcode.wdl" as Split
import "MergeSortBam.wdl" as Merge
import "CreateCountMatrix.wdl" as Count
import "StarAlignBamSingleEnd.wdl" as StarAlignBam
import "TagGeneExon.wdl" as TagGeneExon
import "SequenceDataWithMoleculeTagMetrics.wdl" as Metrics
import "TagSortBam.wdl" as TagSortBam
import "RunEmptyDrops.wdl" as RunEmptyDrops
import "ZarrUtils.wdl" as ZarrUtils
import "Picard.wdl" as Picard
import "MarkDuplicates.wdl" as MarkDuplicates

workflow Optimus {
  meta {
    description: "The optimus 3' pipeline processes 10x genomics sequencing data based on the v2 chemistry. It corrects cell barcodes and UMIs, aligns reads, marks duplicates, and returns data as alignments in BAM format and as counts in sparse matrix exchange format."
  }
  # version of this pipeline
  String version = "optimus_v0.4.0"

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

  # whether to convert the outputs to Zarr format, by default it's set to true
  Boolean output_zarr = true

  Int preemptible = 3

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
    output_zarr: "whether to run the taks that converts the outputs to Zarr format, by default it's true"
    preemptible: "(optional) if non-zero, request a pre-emptible instance and allow for this number of preemptions before running the task on a non preemptible machine"
  }

  scatter (index in indices) {
    call FastqToUBam.FastqToUBam {
      input:
        fastq_file = r2_fastq[index],
        sample_id = sample_id,
        fastq_suffix = fastq_suffix,
        preemptible = preemptible
    }

    # if the index is passed, attach it to the bam file
    if (defined(i1_fastq)) {
      Array[File] non_optional_i1_fastq = select_first([i1_fastq])
      call Attach.Attach10xBarcodes as AttachBarcodes {
        input:
          r1_fastq = r1_fastq[index],
          i1_fastq = non_optional_i1_fastq[index],
          r2_unmapped_bam = FastqToUBam.bam_output,
          whitelist = whitelist,
          preemptible = preemptible
      }
    }

    # if the index is not passed, proceed without it.
    if (!defined(i1_fastq)) {
      call Attach.Attach10xBarcodes as AttachBarcodesNoIndex {
        input:
          r1_fastq = r1_fastq[index],
          r2_unmapped_bam = FastqToUBam.bam_output,
          whitelist = whitelist,
          preemptible = preemptible
      }
    }

    File barcoded_bam = select_first([AttachBarcodes.bam_output, AttachBarcodesNoIndex.bam_output])
  }

  call Merge.MergeSortBamFiles as MergeUnsorted {
    input:
      bam_inputs = barcoded_bam,
      sort_order = "unsorted",
      preemptible = preemptible
  }

  call Split.SplitBamByCellBarcode {
    input:
      bam_input = MergeUnsorted.output_bam,
      preemptible = preemptible
  }

  scatter (bam in SplitBamByCellBarcode.bam_output_array) {
    call StarAlignBam.StarAlignBamSingleEnd as StarAlign {
      input:
        bam_input = bam,
        tar_star_reference = tar_star_reference,
        preemptible = preemptible
    }

    call TagGeneExon.TagGeneExon as TagGenes {
      input:
        bam_input = StarAlign.bam_output,
        annotations_gtf = annotations_gtf,
        preemptible = preemptible
    }

    call Picard.SortBamAndIndex as PreUMISort {
      input:
        bam_input = TagGenes.bam_output,
        preemptible = preemptible
    }

    call MarkDuplicates.MarkDuplicatesUmiTools as MarkDuplicates {
      input:
        bam_input = PreUMISort.bam_output,
        bam_index = PreUMISort.bam_index,
        preemptible = preemptible
    }

    call Picard.SortBamAndIndex as PostUMISort {
      input:
        bam_input = MarkDuplicates.bam_output,
        preemptible = preemptible
    }

    call TagSortBam.GeneSortBam {
      input:
        bam_input = MarkDuplicates.bam_output,
        preemptible = preemptible
    }

    call TagSortBam.CellSortBam {
      input:
        bam_input = MarkDuplicates.bam_output,
        preemptible = preemptible
    }

    call Metrics.CalculateGeneMetrics {
      input:
        bam_input = GeneSortBam.bam_output,
        preemptible = preemptible
    }

    call Metrics.CalculateCellMetrics {
      input:
        bam_input = CellSortBam.bam_output,
        preemptible = preemptible
    }

    call Picard.SortBam as PreCountSort {
      input:
        bam_input = MarkDuplicates.bam_output,
        sort_order = "queryname",
        preemptible = preemptible
    }

    call Count.CreateSparseCountMatrix {
      input:
        bam_input = PreCountSort.bam_output,
        gtf_file = annotations_gtf,
        preemptible = preemptible
    }
  }

  call Merge.MergeSortBamFiles as MergeSorted {
    input:
      bam_inputs = PostUMISort.bam_output,
      sort_order = "coordinate",
      preemptible = preemptible
  }

  call Metrics.MergeGeneMetrics {
    input:
      metric_files = CalculateGeneMetrics.gene_metrics,
      preemptible = preemptible
  }

  call Metrics.MergeCellMetrics {
    input:
      metric_files = CalculateCellMetrics.cell_metrics,
      preemptible = preemptible
  }

  call Count.MergeCountFiles {
    input:
      sparse_count_matrices = CreateSparseCountMatrix.sparse_count_matrix,
      row_indices = CreateSparseCountMatrix.row_index,
      col_indices = CreateSparseCountMatrix.col_index,
      preemptible = preemptible
  }

  call RunEmptyDrops.RunEmptyDrops {
    input:
      sparse_count_matrix = MergeCountFiles.sparse_count_matrix,
      row_index = MergeCountFiles.row_index,
      col_index = MergeCountFiles.col_index,
      preemptible = preemptible
  }

  if (output_zarr) {
    call ZarrUtils.OptimusZarrConversion{
      input:
        sample_id=sample_id,
        cell_metrics = MergeCellMetrics.cell_metrics,
        gene_metrics = MergeGeneMetrics.gene_metrics,
        sparse_count_matrix = MergeCountFiles.sparse_count_matrix,
        cell_id = MergeCountFiles.row_index,
        gene_id = MergeCountFiles.col_index,
        empty_drops_result = RunEmptyDrops.empty_drops_result,
        preemptible = preemptible
    }
}

  output {
      # version of this pipeline
      String pipeline_version = version

      File bam = MergeSorted.output_bam
      File matrix = MergeCountFiles.sparse_count_matrix
      File matrix_row_index = MergeCountFiles.row_index
      File matrix_col_index = MergeCountFiles.col_index
      File cell_metrics = MergeCellMetrics.cell_metrics
      File gene_metrics = MergeGeneMetrics.gene_metrics
      File cell_calls = RunEmptyDrops.empty_drops_result

      # zarr
      Array[File]? zarr_output_files = OptimusZarrConversion.zarr_output_files
  }
}
