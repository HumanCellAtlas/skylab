import "BuildStarReferenceBundle.wdl" as buildStarRef

workflow buildStarRefs {
  ## meta file list the following info
  ## name1: referencec name
  ## name2: reference alias name
  ## fasta: reference fastafile name
  ## gtf: annotaiton reference file
  File ref_metadata_file
  Array[Object] metadata = read_objects(ref_metadata_file)
  
  scatter (record in metadata){
    String ref_name = record.name1
    String name2 = record.name2
    File fasta = record.fasta
    File gtf = record.gtf
    
    call buildStarRef.BuildStarReference as StarRef {
      input:
        ref_fasta = fasta,
        gtf_file = gtf,
        ref_name = ref_name
    }
  }
}
