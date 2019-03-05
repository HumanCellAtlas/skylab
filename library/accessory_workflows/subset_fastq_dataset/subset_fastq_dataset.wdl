task SubsetFastqDatasetTask {
     String star_align_bucket
     String fastq_bucket
     String keep_region

     command {
       mkdir output/
       prepSubset.sh -a ${star_align_bucket} -f ${fastq_bucket} -r ${keep_region} -o output/
       tar cvzf output.tar.gz output/
     }

     runtime {
       docker: "quay.io/humancellatlas/secondary-analysis-subset-fastq:0.0.1"
       memory: "28 GB"
       disks: "local-disk 3000 HDD"
       cpu: "8"
     }

     output {
       File output_tar = "output.tar.gz"
     }
}

workflow SubsetFastqDataset {
     String alignment_location
     String fastq_location
     String keepregion

     call SubsetFastqDatasetTask {
         input:
	    star_align_bucket = alignment_location,
	    fastq_bucket = fastq_location,
	    keep_region = keepregion
     }

     output {
       File output_tar = SubsetFastqDatasetTask.output_tar
     }
}