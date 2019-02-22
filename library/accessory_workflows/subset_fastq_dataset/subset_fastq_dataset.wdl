task SubsetFastqDataset {
     String star_align_bucket
     String fastq_bucket
     String keep_region

     command {
       mkdir output/
       prepSubset.sh -a ${star_align_bucket} -f ${fastq_bucket} -r ${keep_region} -o output
       tar cvzf output.tar.gz output/
     }

     runtime {
       docker: ""
       memory: "28 GB"
       disks: "local-disk 100 HDD"
       cpu: "8"
     }

     output {
       File output_tar = "output.tar.gz"
     }
}
