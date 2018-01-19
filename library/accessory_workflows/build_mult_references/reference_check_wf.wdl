import "star.wdl" as RunStar
import "featurecounts.wdl" as fc

workflow CheckReference {
  ## meta file has list info
  ## sraID: sraID for input file
  ## star_genome: star genome name, such GRCh38_RefSeq, CRCh38_GencodeV27
  ## gtf: annotation file names
  ## length: fastq reads length, 25 or 100
  File input_metadata
  Array[Object] metadata = read_objects(input_metadata)
  String sra_dir

  scatter(record in metadata) {
    call RunStar.Star as Star {  
      input:
        input_fastq_read1 = sra_dir+record.sraID+"_1.fastq.gz",
        input_fastq_read2 = sra_dir+record.sraID+"_2.fastq.gz",
        gtf = record.gtf,
        star_genome = record.star_genome,
        sample_tag = record.sraID+"_"+record.length+"_"+record.ref_name,
        pu_tag = record.sraID+"_"+record.length+"_"+record.ref_name,
        lib_tag = record.sraID+"_"+record.length+"_"+record.ref_name,
        id_tag = record.sraID+"_"+record.length+"_"+record.ref_name
    }

    call fc.FeatureCountsUniqueMapping as uniqcount {
      input:
        aligned_bam = Star.output_bam,
        gtf = record.gtf,
        fc_out = record.sraID+"_"+record.length+"_"+record.ref_name
    }

    call fc.FeatureCountsMultiMapping as multcount {
      input:
        aligned_bam=Star.output_bam,
        gtf = record.gtf,
        fc_out = record.sraID+"_"+record.length+"_"+record.ref_name
    }
 }

 output {
   Star.logs
   Star.junction_table
   Star.output_bam
   Star.output_bam_trans
   uniqcount.genes
   uniqcount.exons
   uniqcount.trans
   multcount.genes
   multcount.exons
   multcount.trans
  }
}
