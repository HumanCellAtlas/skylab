import "FastqToUBam.wdl" as fq2bam
import "Attach10xBarcodes.wdl" as attach
import "SplitBamByCellBarcode.wdl" as split
import "CollectMultiplePicardMetrics.wdl" as collect
import "MergeSortBam.wdl" as merge
import "CreateCountMatrix.wdl" as count
import "AlignTagCorrectUmis.wdl" as AlignTagCorrectUmis

# The optimus 3' pipeline processes 10x genomics sequencing data based on the v2
# chemistry. It corrects cell barcodes and UMIs, aligns reads, marks duplicates, and
# returns data as alignments in BAM format and as counts in sparse matrix exchange
# format.
workflow Optimus {

  Array[Array[File]] fastq_inputs  # array of arrays of matched fastqs [ [r1, r2, i1], ... ]
  File whitelist  # 10x genomics cell barcode whitelist for 10x V2
  File tar_star_reference  # star reference
  File annotations_gtf  # gtf containing annotations for gene tagging
  File ref_genome_fasta  # genome fasta file
  String sample_id  # name of sample matching this file, inserted into read group header
  String fastq_suffix = ""  # when running in green box, need to add ".gz" for picard to detect

  # this scatters matched [r1, r2, i1] fastq arrays
  scatter (fastqs in fastq_inputs) {
    call fq2bam.FastqToUBam {
      input:
        fastq_file = fastqs[1],
        sample_id = sample_id,
        fastq_suffix = fastq_suffix
    }

    call attach.Attach10xBarcodes {
      input:
        r1 = fastqs[0],
        i1 = fastqs[2],
        u2 = FastqToUBam.bam_output,
        whitelist = whitelist
    }

    call split.SplitBamByCellBarcode {
      input:
        bam_input = Attach10xBarcodes.bam_output
    }

    call AlignTagCorrectUmis.AlignTagCorrectUmis {
      input:
        bam_array = SplitBamByCellBarcode.bam_output_array,
        tar_star_reference = tar_star_reference,
        annotations_gtf = annotations_gtf
    }
  }

  call merge.MergeSortBamFiles {
    input:
      bam_inputs = AlignTagCorrectUmis.bam_outputs
  }

  call count.DropSeqToolsDigitalExpression {
    input:
      bam_input = MergeSortBamFiles.output_bam,
      whitelist = whitelist
  }

  call collect.CollectMultipleMetrics {
    input:
      aligned_bam = MergeSortBamFiles.output_bam,
      ref_genome_fasta = ref_genome_fasta,
      output_filename = sample_id
  }

  output {
      File bam = MergeSortBamFiles.output_bam
      File matrix = DropSeqToolsDigitalExpression.matrix_output
      File matrix_summary = DropSeqToolsDigitalExpression.matrix_summary
      Array[Array[File]] tag_gene_exon_log = AlignTagCorrectUmis.tag_gene_exon_log
      Array[Array[File]] umi_metrics = AlignTagCorrectUmis.umi_metrics
      Array[Array[File]] duplicate_metrics = AlignTagCorrectUmis.duplicate_metrics
      File picard_metrics = CollectMultipleMetrics.alignment_metrics
  }
}
