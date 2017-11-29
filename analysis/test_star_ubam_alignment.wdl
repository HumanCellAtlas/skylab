import "StarAlignFastqSingleEnd.wdl" as star_fastq
import "StarAlignBamSingleEnd.wdl" as star_bam
import "FastqToUBam.wdl" as fq2bam
import "Attach10xBarcodes.wdl" as attach

workflow test_ubam_alignment {
  File r1  # forward fastq file; contains umi, cell barcode
  File r2  # reverse fastq file; contains alignable genomic information
  File i1  # index fastq file; contains sample barcode
  String sample_name  # name of sample matching this file, inserted into read group header
  File tar_star_reference  # tarball containing star reference.

  call fq2bam.FastqToUBam {
    input:
      gz_fastq_file = r2,
      sample_name = sample_name
  }

  call attach.Attach10xBarcodes {
    input:
      r1 = r1,
      i7 = i1,
      u2 = FastqToUBam.bam
  }

  call star_bam.StarAlignBamSingleEnd {
    input:
      bam_input = Attach10xBarcodes.barcoded_bam,
      tar_star_reference = tar_star_reference
  }

  call star_fastq.StarAlignFastqSingleEnd {
    input:
      fastq_input = r2,
      tar_star_reference = tar_star_reference
  }
  
  output {
    File bam_to_bam = StarAlignBamSingleEnd.bam
    File bam_alignment_log = StarAlignBamSingleEnd.alignment_log
    File fastq_to_bam = StarAlignFastqSingleEnd.bam
    File fastq_alignment_log = StarAlignFastqSingleEnd.alignment_log
  }
}
