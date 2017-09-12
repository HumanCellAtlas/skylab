# stage SETUP_CHUNKS(
#     in  string sample_id,
#     in  map[]  sample_def,
#     in  string chemistry_name,
#     in  map    custom_chemistry_def,
#     out map[]  chunks,
#     out map    chemistry_def,
#     out string barcode_whitelist,
#     src py     "stages/common/setup_chunks",
# )
task setup_chunks {

  String sample_id
  File sample_def
  Array[File] read_paths
  
  String chemistry_name = "auto"
  String d = "$"
  String l = "{"
  String r = "}"

  command {

    # The sample_def maps contain a key "read_path" that refers to a path where SETUP_CHUNKS is
    # supposed to find the input fastq files. If that path is hidden within the Map, cromwell
    # will never know about it, so we have to expose it via a separate input, read_path. But,
    # the path in the Map is going to be different from the path in the Array[File] because the
    # latter is going to have a bunch of cromwell stuff prepended to it. So we need to put those
    # in a location that SETUP_CHUNKS is expecting.
    set -x
    for read_path in '${sep="\' \'" read_paths}'; do
      mkdir -p read_paths/${d}${l}read_path##*call-setup_chunks/inputs${r}
      cp -al "${d}read_path" read_paths/${d}${l}read_path##*call-setup_chunks/inputs${r}
    done
    
    sample_def_string=${d}(<${sample_def})
    /opt/tenx/run_in_10x_env.bash martian stage run SETUP_CHUNKS main \
      --sample_id '${sample_id}' \
      --chemistry_name '${chemistry_name}' \
      --sample_def "${d}sample_def_string"
  }
  

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File chunks = "chunks" 
    File chemistry_def = "chemistry_def"
    String barcode_whitelist = read_string("barcode_whitelist")
  }
}

# stage CHUNK_READS(
#     in  map[] chunks,
#     in  int   reads_per_file,
#     out map[] out_chunks,
#     src py    "stages/common/chunk_reads",
# ) split using (
#     in  map   read_chunks,
# )
task chunk_reads_split {
  File chunks
  Int reads_per_file
  String d = "$"
  String l = "{"
  String r = "}"

  command {
    set -x
    
    # chunks is a file with json. Read it into a string before passing it to Martian.
    chunks_string=${d}(<${chunks})

    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS split \
      --chunks "${d}chunks_string" \
      --reads_per_file '${reads_per_file}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] split_files = glob("martian_split_*")
  }
}

task chunk_reads_main {
  File chunks 
  Int reads_per_file
  File split_file
  Array[File] read_paths

  String d = "$"
  String l = "{"
  String r = "}"

  command {
    set -x
    
    # Move reads paths to $PWD/read_paths so string refer to actual paths
    for read_path in '${sep="\' \'" read_paths}'; do
      mkdir -p read_paths/${d}${l}read_path##*call-chunk_reads_main/shard-*/inputs${r}
      cp -al "${d}read_path" read_paths/${d}${l}read_path##*call-chunk_reads_main/shard-*/inputs${r}
    done


    chunks_string=${d}(<${chunks})

    # This stage is going to call an executable that writes to a "files" directory    
    mkdir files

    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS main \
      --chunks "${d}chunks_string" \
      --reads_per_file '${reads_per_file}' \
      --split_file '${split_file}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File out_chunks = "out_chunks"
    Array[File] fastq_chunks = glob("files/*.fastq")
  }
}

task chunk_reads_join {
  Array[File] out_chunks_output

  String d = "$"
  String l = "{"
  String r = "}"

  command {
    set -x

    # Merge the array of chunks jsons written by each main step into a single file,
    # and read that to a string
    jq -s '.' '${sep="\' \'" out_chunks_output}' > merged_out_chunks
    merged_out_chunks_string=${d}(<merged_out_chunks)
    
    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS join \
      --out_chunks_output "${d}merged_out_chunks_string"
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File out_chunks = "out_chunks"
  }
}

# stage EXTRACT_READS(
#     in  map[]    chunks,
#     in  map      chemistry_def,
#     in  string   barcode_whitelist,
#     in  int      reads_per_file,
#     in  float    subsample_rate,
#     in  int      initial_reads,
#     in  map[]    primers,
#     in  map      align,
#     in  int      rna_read_length,
#     in  bool     skip_metrics,
#     out pickle   chunked_reporter,
#     out json     summary,
#     out json     barcode_counts,
#     out fastq[]  reads,
#     out fastq[]  read2s,
#     out bam[]    trimmed_seqs,
#     out int[]    gem_groups,
#     out string[] read_groups,
#     out map      align,
#     out string[] bam_comments,
#     src py       "stages/common/extract_reads",
# ) split using (
#     in  map      read_chunks,
#     in  int      gem_group,
#     in  bool     reads_interleaved,
#     in  bool     barcode_rc,
# )
task extract_reads_split {
  File chunks
  String barcode_whitelist
  
  String d = "$"

  command {

    chunks_string=${d}(<${chunks})

    /opt/tenx/run_in_10x_env.bash martian stage run EXTRACT_READS split \
      --chunks "${d}chunks_string" \
      --barcode_whitelist ${barcode_whitelist}
  }
 
  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] split_files = glob("martian_split_*")
  }
}

task extract_reads_main {
  File chunks
  File chemistry_def
  String barcode_whitelist
  Int reads_per_file
  Float subsample_rate
  Array[Map[String, String]] primers

  File split_file
  Array[File] fastq_chunks
 
  String d = '$' 
  String l = '{'
  String r = '}'
 
  command {

    # Put the fastq file into the $PWD/files directory where martian is expecting
    # them
    mkdir files
    for fastq_chunk in ${sep=" " fastq_chunks}; do
      ln "${d}fastq_chunk" files/"${d}(basename ${d}fastq_chunk)"
    done
    
    # Read the chunks and chemistry json files to strings so we can pass them as
    # a Martian map
    chunks_string=${d}(<${chunks})
    chemistry_def_string=${d}(<${chemistry_def})

    /opt/tenx/run_in_10x_env.bash martian stage run EXTRACT_READS main \
      --chunks "${d}chunks_string" \
      --chemistry_def "${d}chemistry_def_string" \
      --barcode_whitelist ${barcode_whitelist} \
      --reads_per_file '${reads_per_file}' \
      --subsample_rate '${subsample_rate}' \
      --skip_metrics no \
      --primers '[${sep="," primers}]' \
      --split_file '${split_file}'
    
    # The Martian cli adapter dumps these outputs as json file, so convert them
    # to a format Cromwell can parse
    jq -r '.[]' read_groups > read_groups.lines
    jq -r '.[]' gem_groups > gem_groups.lines
    jq -r '.[]' bam_comments > bam_comments.lines

  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] trimmed_seqs = glob("trimmed_seqs.bam/*")
    Array[File] reads = glob("reads.fastq/*")
    File barcode_counts = "barcode_counts.json"
    File chunked_reporter = "chunked_reporter.pickle"
    Array[String] read_groups = read_lines("read_groups.lines")
    Array[Int] gem_groups = read_lines("gem_groups.lines")
    Array[String] bam_comments = read_lines("bam_comments.lines")
    File align = "align"
  }
}

task extract_reads_join {
  String align 
  String barcode_whitelist
  File chemistry_def
  Array[Array[File]] reads_output
  Array[Array[File]] trimmed_seqs_output
  Array[Array[String]] bam_comments_output
  Array[Array[Int]] gem_groups_output
  Array[Array[String]] read_groups_output
  Array[File] chunked_reporter_output
  Array[File] barcode_counts_output
    
  Array[File] split_files
  
  String d = "$"
  String l = "{"
  String r = "}"

  command {
    set -x  
    
    # There's a lot here, but it's mostly the same thing repeated. All the input reads and trimmed_seq
    # files need to be moved to paths in $PWD. Also, we have to prepare a json list of lists that
    # has all these new paths so we can pass it to martian. The jq ' |= .+ [ ]' syntax just appends 
    # something to a list. Finally, we need to create a list of empty lists called read2s that's the
    # same length as the list of reads.
    idx=0
    reads_arr_arg=[]
    read2s_arr_arg=[]
    for reads_arr in ${d}(echo '[${sep="," reads_output}]' | jq -c '.[]'); do
      mkdir -p out_reads/${d}idx
      reads_arg=[]
      for reads_path in ${d}(echo ${d}reads_arr | jq -r '.[]'); do
        new_path=out_reads/${d}idx/${d}(basename ${d}reads_path)
        ln ${d}reads_path ${d}new_path
        reads_arg=${d}( echo ${d}reads_arg | jq -c --arg v ${d}new_path '. |= .+ [$v]')
      done
      reads_arr_arg=${d}( echo ${d}reads_arr_arg | jq -c --arg v ${d}reads_arg '. |= .+ [$v | fromjson]')
      # Create a dummy read2s_output arg    
      read2s_arg=[]
      read2s_arr_arg=${d}( echo ${d}read2s_arr_arg | jq -c --arg v ${d}read2s_arg '. |= .+ [$v | fromjson]')
      idx=${d}((idx+1))
    done

    
    idx=0
    trimmed_seqs_arr_arg=[]
    for trimmed_seqs_arr in ${d}(echo '[${sep="," trimmed_seqs_output}]' | jq -c '.[]'); do
      mkdir -p out_trimmed_seqs/${d}idx
      trimmed_seqs_arg=[]
      for trimmed_seqs_path in ${d}(echo ${d}trimmed_seqs_arr | jq -r '.[]'); do
        new_path=out_trimmed_seqs/${d}idx/${d}(basename ${d}trimmed_seqs_path)
        ln ${d}trimmed_seqs_path ${d}new_path
        trimmed_seqs_arg=${d}( echo ${d}trimmed_seqs_arg | jq -c --arg v ${d}new_path '. |= .+ [$v]')
      done
      trimmed_seqs_arr_arg=${d}( echo ${d}trimmed_seqs_arr_arg | jq -c --arg v ${d}trimmed_seqs_arg '. |= .+ [$v | fromjson]')
      idx=${d}((idx+1))
    done
  
    # The gem_groups we want are stored in each of the files produced by the split step, so pull them out
    # and stick them into a big list of lists of ints.
    gem_groups=["${d}(jq -s '.[].gem_group | tostring' ${sep=' ' split_files} | jq -rs 'join(",")')"]
    chemistry_def_string=${d}(<${chemistry_def})
    
    /opt/tenx/run_in_10x_env.bash martian stage run EXTRACT_READS join \
      --align '${align}' \
      --barcode_whitelist ${barcode_whitelist} \
      --chemistry_def "${d}chemistry_def_string" \
      --gem_group_split "${d}gem_groups" \
      --reads_output "${d}reads_arr_arg" \
      --read2s_output "${d}read2s_arr_arg" \
      --trimmed_seqs_output "${d}trimmed_seqs_arr_arg" \
      --bam_comments_output '[${sep=',' bam_comments_output}]' \
      --gem_groups_output '[${sep=',' gem_groups_output}]' \
      --read_groups_output '[${sep=',' read_groups_output}]' \
      --chunked_reporter_output '["${sep='","' chunked_reporter_output}"]' \
      --barcode_counts_output '["${sep='","' barcode_counts_output}"]'
    
    # As above, convert json to Cromwell formats
    jq -r '.[]' read_groups > read_groups.lines
    jq -r '.[]' gem_groups > gem_groups.lines
    jq -r '.[]' bam_comments > bam_comments.lines
    jq -r '.[]' reads > reads.lines
    jq -r '.[]' trimmed_seqs > trimmed_seqs.lines
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File summary = "summary.json"
    File barcode_counts = "barcode_counts.json"
    Array[String] read_groups = read_lines("read_groups.lines")
    Array[Int] gem_groups = read_lines("gem_groups.lines")
    Array[String] bam_comments = read_lines("bam_comments.lines")
    Array[String] reads = read_lines("reads.lines")
    Array[String] trimmed_seqs = read_lines("trimmed_seqs.lines")
    File out_align = "align"
  }
}

task align_reads_split {
  Array[String] read1s
  Array[String] read_groups
  File reference_path
    
  String d = "$"

  command {
    /opt/tenx/run_in_10x_env.bash martian stage run ALIGN_READS split \
      --reads '["${sep='","' read1s}"]' \
      --read2s '[]' \
      --read_groups '["${sep='","' read_groups}"]' \
      --reference_path '${reference_path}' \
      --threads ${d}(nproc)
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] split_files = glob("martian_split_*")
  }
}

task align_reads_main {
  Int max_hits_per_read
  File reference_path
  Array[File] read_chunk
  File split_file
    
  String d = "$"
  
  # Martian passes these around as arrays but the code seems to expect that it's
  # just a single file. And it looks like the array has just one element, so just
  # pull that out here.
  File first_read_chunk = select_first(read_chunk)

  command {
    set -x
    read_group=${d}(jq -r '.read_group' ${split_file})

    /opt/tenx/run_in_10x_env.bash martian stage run ALIGN_READS main \
      --threads ${d}(nproc) \
      --read_chunk '${first_read_chunk}' \
      --reference_path '${reference_path}' \
      --max_hits_per_read ${max_hits_per_read} \
      --read_group "${d}read_group"
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File genome_output = "genome_output.bam"
  }
}


task attach_bcs_and_umis_main {
  File reference_path
  String barcode_whitelist
  File barcode_counts
  File chemistry_def
  Array[Int] gem_groups
  Int umi_min_qual_threshold
  Array[String] bam_comments
  
  # See the notes in workflow below, but here we're just unpacking a Pair[X, Pair[Y, Z]]
  # created by a nested zip
  Pair[Int, Pair[Array[File], File]] chunk
  Int gem_group = chunk.left
  Pair[Array[File], File] chunk_right = chunk.right
  Array[File] trimmed_seqs_arr = chunk_right.left
  File chunk_trimmed_input = select_first(trimmed_seqs_arr)
  File chunk_genome_input = chunk_right.right

  String d = "$"

  command {
    
    chemistry_def_string=${d}(<${chemistry_def})   

    /opt/tenx/run_in_10x_env.bash martian stage run ATTACH_BCS_AND_UMIS main \
     --reference_path '${reference_path}' \
     --barcode_whitelist ${barcode_whitelist} \
     --barcode_counts '${barcode_counts}' \
     --gem_groups '[${sep=',' gem_groups}]' \
     --gem_group '${gem_group}' \
     --umi_min_qual_threshold '${umi_min_qual_threshold}' \
     --barcode_counts '${barcode_counts}' \
     --chunk_genome_input '${chunk_genome_input}' \
     --chunk_trimmed_input '${chunk_trimmed_input}' \
     --paired_end false \
     --skip_metrics false \
     --correct_barcode true \
     --rescue_multimappers true \
     --barcode_confidence_threshold 0.975 \
     --bam_comments '["${sep='","' bam_comments}"]' \
     --chemistry_def "${d}chemistry_def_string"

  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File bam_output = "output.bam"
    Int num_alignments = read_int("num_alignments")
    File chunked_reporter = "chunked_reporter.pickle"
  }
}

task attach_bcs_and_umis_join {
  File reference_path
  Array[Int] num_alignments_output
  Array[File] output_output
  Array[File] chunked_reporter_output

  command {
    /opt/tenx/run_in_10x_env.bash martian stage run ATTACH_BCS_AND_UMIS join \
     --reference_path '${reference_path}' \
     --num_alignments_output '[${sep=',' num_alignments_output}]' \
     --output_output '["${sep='","' output_output}"]' \
     --chunked_reporter_output '["${sep='","' chunked_reporter_output}"]'

    jq -r '.[]' num_alignments > num_alignments.lines
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File summary = "summary.json"
    File barcode_summary = "barcode_summary.h5"
    Array[Int] num_alignments = read_lines("num_alignments.lines")
  }
}

task bucket_by_bc_split {
  Array[File] inputs
  Array[Int] num_alignments

  command {
    /opt/tenx/run_in_10x_env.bash martian stage run BUCKET_BY_BC split \
     --num_alignments '[${sep=',' num_alignments}]' \
     --inputs '["${sep='","' inputs}"]'

    jq -c '.read_groups' martian_split_0 > read_groups
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File read_groups = "read_groups"
  }
}

task bucket_by_bc_main {
  File chunk_input
  File read_groups 
  
  String d = "$"

  command {
    read_groups_string=${d}(<${read_groups})

    mkdir files

    /opt/tenx/run_in_10x_env.bash martian stage run BUCKET_BY_BC main \
      --chunk_input '${chunk_input}' \
      --read_groups "${d}read_groups_string" \
      --nbases 2
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }
  
  # This step splits everything up by pairs of nucleotides. Each pair has to be handled separately,
  # so just break them out verbosely here. Hopefully that's clearer than trying to be clever about
  # it.
  output {
    File buckets = "buckets"
    File AA_bam = "files/bc_AA.bam"
    File AT_bam = "files/bc_AT.bam"
    File AG_bam = "files/bc_AG.bam"
    File AC_bam = "files/bc_AC.bam"
    File TA_bam = "files/bc_TA.bam"
    File TT_bam = "files/bc_TT.bam"
    File TG_bam = "files/bc_TG.bam"
    File TC_bam = "files/bc_TC.bam"
    File GA_bam = "files/bc_GA.bam"
    File GT_bam = "files/bc_GT.bam"
    File GG_bam = "files/bc_GG.bam"
    File GC_bam = "files/bc_GC.bam"
    File CA_bam = "files/bc_CA.bam"
    File CT_bam = "files/bc_CT.bam"
    File CG_bam = "files/bc_CG.bam"
    File CC_bam = "files/bc_CC.bam"
    File null_bam = "files/bc_.bam"
  }
}

task sort_by_bc_main {
  Array[File] bucketed_bams
  
  String d = "$"

  command {
    set -x
    
    # Like above, move the bam inputs to PWD and prepare a json list of the new paths
    idx=0
    bucket=[]
    for bam_path in '${sep="\' \'" bucketed_bams}'; do
      mkdir -p bucketed_bams/${d}idx
      new_path=bucketed_bams/$idx/"${d}(basename ${d}bam_path)"
      ln "${d}bam_path" ${d}new_path
      bucket=${d}(echo ${d}bucket | jq -c --arg v ${d}new_path '. |= .+ [$v]')
      idx=${d}((idx+1))
    done
    
    mkdir files

    /opt/tenx/run_in_10x_env.bash martian stage run SORT_BY_BC main \
      --bucket "${d}bucket"
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File default = "default.bam"
    Int total_reads = read_int("total_reads")
  }
}

task sort_by_bc_join {
  File sorted_AA
  File sorted_AT
  File sorted_AC
  File sorted_AG
  File sorted_TA
  File sorted_TT
  File sorted_TC
  File sorted_TG
  File sorted_CA
  File sorted_CT
  File sorted_CC
  File sorted_CG
  File sorted_GA
  File sorted_GT
  File sorted_GC
  File sorted_GG
  File sorted_null
  Int total_AA
  Int total_AT
  Int total_AC
  Int total_AG
  Int total_TA
  Int total_TT
  Int total_TC
  Int total_TG
  Int total_CA
  Int total_CT
  Int total_CC
  Int total_CG
  Int total_GA
  Int total_GT
  Int total_GC
  Int total_GG
  Int total_null

  command {
    set -x
    
    mkdir files
    
    # Yikes
    /opt/tenx/run_in_10x_env.bash martian stage run SORT_BY_BC join \
      --prefix_split '["AA", "AT", "AC", "AG", "TA", "TT", "TC", "TG", "CA", "CT", "CC", "CG", "GA", "GT", "GC", "GG", ""]' \
      --default_output '["${sorted_AA}", "${sorted_AT}", "${sorted_AC}", "${sorted_AG}","${sorted_TA}", "${sorted_TT}", "${sorted_TC}", "${sorted_TG}","${sorted_CA}", "${sorted_CT}", "${sorted_CC}", "${sorted_CG}","${sorted_GA}", "${sorted_GT}", "${sorted_GC}", "${sorted_GG}", "${sorted_null}"]' \
      --total_reads_output '[${total_AA}, ${total_AT}, ${total_AC}, ${total_AG},${total_TA}, ${total_TT}, ${total_TC}, ${total_TG},${total_CA}, ${total_CT}, ${total_CC}, ${total_CG},${total_GA}, ${total_GT}, ${total_GC}, ${total_GG}, ${total_null}]'
      
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File default = "default.bam"
    Int total_reads = read_int("total_reads")
  }
}


task mark_duplicates_split {
  File input_bam
  
  command {
    /opt/tenx/run_in_10x_env.bash martian stage run MARK_DUPLICATES split \
      --input '${input_bam}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] split_files = glob("martian_split_*")
  }

}

task mark_duplicates_main {

  File input_bam
  File align
  File split_file
  
  String d = "$"

  command {

    align_string=${d}(<${align})

    /opt/tenx/run_in_10x_env.bash martian stage run MARK_DUPLICATES main \
      --input '${input_bam}' \
      --split_file '${split_file}' \
      --align "${d}align_string"
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File chunked_reporter = "chunked_reporter.pickle"
    File output_bam = "output.bam"
  }
}

task mark_duplicates_join {
  Array[File] output_output
  Array[File] chunked_reporter_output

  command {
    /opt/tenx/run_in_10x_env.bash martian stage run MARK_DUPLICATES join \
      --output_output '["${sep='","' output_output}"]' \
      --chunked_reporter_output '["${sep='","' chunked_reporter_output}"]'
  }
  
  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File summary = "summary.json"
  }

}

task count_genes_split {
  Array[File] bam_inputs
  String barcode_whitelist
  Array[Int] gem_groups
  File reference_path
  File barcode_summary
  File chemistry_def

  String d = "$" 

  command {

    chemistry_def_string=${d}(<${chemistry_def})   

    idx=0
    linked_bam_inputs=[]
    for bam_path in '${sep="\' \'" bam_inputs}'; do
      mkdir -p bam_inputs/${d}idx
      new_path=bam_inputs/$idx/"${d}(basename ${d}bam_path)"
      ln "${d}bam_path" ${d}new_path
      linked_bam_inputs=${d}(echo ${d}linked_bam_inputs | jq -c --arg v ${d}new_path '. |= .+ [$v]')
      idx=${d}((idx+1))
    done

    /opt/tenx/run_in_10x_env.bash martian stage run COUNT_GENES split \
      --inputs "${d}linked_bam_inputs" \
      --barcode_whitelist ${barcode_whitelist} \
      --gem_groups '[${sep=',' gem_groups}]' \
      --reference_path '${reference_path}' \
      --barcode_summary '${barcode_summary}' \
      --chemistry_def "${d}chemistry_def_string"
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] split_files = glob("martian_split_*")
  }
}


task count_genes_main {
  Array[File] bam_inputs
  String barcode_whitelist
  Array[Int] gem_groups
  File reference_path
  File barcode_summary
  File chemistry_def
  File split_file
  File align

  String d = "$" 

  command {

    chemistry_def_string=${d}(<${chemistry_def})   
    align_string=${d}(<${align})   

    idx=0
    linked_bam_inputs=[]
    for bam_path in '${sep="\' \'" bam_inputs}'; do
      mkdir -p bam_inputs/${d}idx
      new_path=bam_inputs/$idx/"${d}(basename ${d}bam_path)"
      ln "${d}bam_path" ${d}new_path
      linked_bam_inputs=${d}(echo ${d}linked_bam_inputs | jq -c --arg v ${d}new_path '. |= .+ [$v]')
      idx=${d}((idx+1))
    done

    /opt/tenx/run_in_10x_env.bash martian stage run COUNT_GENES main \
      --inputs "${d}linked_bam_inputs" \
      --barcode_whitelist ${barcode_whitelist} \
      --gem_groups '[${sep=',' gem_groups}]' \
      --reference_path '${reference_path}' \
      --barcode_summary '${barcode_summary}' \
      --chemistry_def "${d}chemistry_def_string" \
      --align "${d}align_string" \
      --split_file '${split_file}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File matrices_h5 = "matrices_h5.h5"
    File chunked_reporter = "chunked_reporter.pickle"
  }
}

task count_genes_join {
  Array[Int] gem_groups
  File chemistry_def
  String sample_id
  Array[File] matrices_h5_output
  Array[File] chunked_reporter_output
  
  String d = "$"

  command {
    chemistry_def_string=${d}(<${chemistry_def})   

    /opt/tenx/run_in_10x_env.bash martian stage run COUNT_GENES join \
      --matrices_h5_output '["${sep='","' matrices_h5_output}"]' \
      --chunked_reporter_output '["${sep='","' chunked_reporter_output}"]' \
      --sample_id '${sample_id}' \
      --gem_groups '[${sep=',' gem_groups}]' \
      --chemistry_def "${d}chemistry_def_string"
  }
  
  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File matrices_h5 = "matrices_h5.h5"
    File matrices_mex = "matrices_mex"
    File reporter_summary = "reporter_summary.json"
    File barcode_summary = "barcode_summary.h5"
  }

}

task filter_barcodes_main {

  String sample_id
  File matrices_h5
  File raw_fastq_summary
  File attach_bcs_summary
  String barcode_whitelist
  Array[Int] gem_groups
  File chemistry_def
  File barcode_summary

  String d = "$"

  command {
    chemistry_def_string=${d}(<${chemistry_def})   

    /opt/tenx/run_in_10x_env.bash martian stage run FILTER_BARCODES main \
      --matrices_h5 '${matrices_h5}' \
      --sample_id '${sample_id}' \
      --gem_groups '[${sep=',' gem_groups}]' \
      --chemistry_def "${d}chemistry_def_string" \
      --raw_fastq_summary '${raw_fastq_summary}' \
      --attach_bcs_summary '${attach_bcs_summary}' \
      --barcode_summary '${barcode_summary}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File summary = "summary.json"
    File filtered_barcodes = "filtered_barcodes.csv"
    File filtered_matrices_h5 = "filtered_matrices_h5.h5"
    File filtered_matrices_mex = "filtered_matrices_mex"
  }
}

workflow count {

  File sample_def
  Array[File] read_paths
  Int reads_per_file
  Float subsample_rate
  String align
  Array[Map[String, String]] primers
  File reference_path
  Int umi_min_qual_threshold
  String sample_id

  call setup_chunks {
    input: sample_id=sample_id, sample_def=sample_def, read_paths=read_paths
  }

  call chunk_reads_split {
    input: chunks=setup_chunks.chunks, reads_per_file=reads_per_file
  }

  scatter(split in chunk_reads_split.split_files) {
    call chunk_reads_main {
      input: split_file=split, chunks=setup_chunks.chunks, reads_per_file=reads_per_file, read_paths=read_paths
    }
  }

  call chunk_reads_join {
    input: out_chunks_output=chunk_reads_main.out_chunks
  }

  call extract_reads_split {
    input: chunks=chunk_reads_join.out_chunks, barcode_whitelist=setup_chunks.barcode_whitelist
  }

  scatter(split_pair in zip(extract_reads_split.split_files, chunk_reads_main.fastq_chunks)) {
      call extract_reads_main {
        input: chunks=chunk_reads_join.out_chunks, barcode_whitelist=setup_chunks.barcode_whitelist,
          primers=primers, reads_per_file=reads_per_file, subsample_rate=subsample_rate,
          chemistry_def=setup_chunks.chemistry_def, split_file=split_pair.left,
          fastq_chunks=split_pair.right
    }
  }

  call extract_reads_join {
    input: align=align, barcode_whitelist=setup_chunks.barcode_whitelist, chemistry_def=setup_chunks.chemistry_def,
      reads_output=extract_reads_main.reads, 
      trimmed_seqs_output=extract_reads_main.trimmed_seqs, bam_comments_output=extract_reads_main.bam_comments,
      gem_groups_output=extract_reads_main.gem_groups, read_groups_output=extract_reads_main.read_groups,
      chunked_reporter_output=extract_reads_main.chunked_reporter, barcode_counts_output=extract_reads_main.barcode_counts,
      split_files=extract_reads_split.split_files
  }

  call align_reads_split {
    input: read1s=extract_reads_join.reads, read_groups=extract_reads_join.read_groups,
      reference_path=reference_path
  }

  scatter(split_pair in zip(align_reads_split.split_files, extract_reads_main.reads)) {
    call align_reads_main {
      input: reference_path=reference_path, split_file=split_pair.left, read_chunk=split_pair.right,
        max_hits_per_read=-1
    }
  }


  # Bear with me here...
  # What I want to do is scatter over zip(gem_groups, trimmed_seqs, genome_output), but in WDL you can only zip two arrays
  # to create Pairs. So here I'm just zip(gem_groups, zip(trimmed_seq, genome_output)) and dealing with it in the task
  Array[Pair[Array[File], File]] trimmed_seqs_and_genome_output = zip(extract_reads_main.trimmed_seqs, align_reads_main.genome_output)
  Array[Pair[Int, Pair[Array[File], File]]] attach_chunks = zip(extract_reads_join.gem_groups, trimmed_seqs_and_genome_output)

  scatter(chunk in attach_chunks) {
    call attach_bcs_and_umis_main {
      input: reference_path=reference_path, barcode_whitelist=setup_chunks.barcode_whitelist, chemistry_def=setup_chunks.chemistry_def,
        gem_groups=extract_reads_join.gem_groups, umi_min_qual_threshold=umi_min_qual_threshold, barcode_counts=extract_reads_join.barcode_counts,
        chunk=chunk, bam_comments=extract_reads_join.bam_comments 
    }
  }

  call attach_bcs_and_umis_join {
    input: reference_path=reference_path, num_alignments_output=attach_bcs_and_umis_main.num_alignments, output_output=attach_bcs_and_umis_main.bam_output,
      chunked_reporter_output=attach_bcs_and_umis_main.chunked_reporter
  }
 
  call bucket_by_bc_split {
    input: inputs=attach_bcs_and_umis_main.bam_output, num_alignments=attach_bcs_and_umis_main.num_alignments
  }

  scatter(chunk_input in attach_bcs_and_umis_main.bam_output) {
    call bucket_by_bc_main {
      input: chunk_input=chunk_input, read_groups=bucket_by_bc_split.read_groups
    }
  }

  call sort_by_bc_main as sort_AA { input: bucketed_bams=bucket_by_bc_main.AA_bam }
  call sort_by_bc_main as sort_AT { input: bucketed_bams=bucket_by_bc_main.AT_bam }
  call sort_by_bc_main as sort_AG { input: bucketed_bams=bucket_by_bc_main.AG_bam }
  call sort_by_bc_main as sort_AC { input: bucketed_bams=bucket_by_bc_main.AC_bam }
  call sort_by_bc_main as sort_TA { input: bucketed_bams=bucket_by_bc_main.TA_bam }
  call sort_by_bc_main as sort_TT { input: bucketed_bams=bucket_by_bc_main.TT_bam }
  call sort_by_bc_main as sort_TG { input: bucketed_bams=bucket_by_bc_main.TG_bam }
  call sort_by_bc_main as sort_TC { input: bucketed_bams=bucket_by_bc_main.TC_bam }
  call sort_by_bc_main as sort_GT { input: bucketed_bams=bucket_by_bc_main.GT_bam }
  call sort_by_bc_main as sort_GG { input: bucketed_bams=bucket_by_bc_main.GG_bam }
  call sort_by_bc_main as sort_GC { input: bucketed_bams=bucket_by_bc_main.GC_bam }
  call sort_by_bc_main as sort_CA { input: bucketed_bams=bucket_by_bc_main.CA_bam }
  call sort_by_bc_main as sort_CT { input: bucketed_bams=bucket_by_bc_main.CT_bam }
  call sort_by_bc_main as sort_CG { input: bucketed_bams=bucket_by_bc_main.CG_bam }
  call sort_by_bc_main as sort_CC { input: bucketed_bams=bucket_by_bc_main.CC_bam }
  call sort_by_bc_main as sort_null { input: bucketed_bams=bucket_by_bc_main.null_bam }
  
  call sort_by_bc_join {
    input:
      sorted_AA=sort_AA.default, sorted_AT=sort_AT.default, sorted_AC=sort_AC.default, sorted_AG=sort_AG.default,
      sorted_TA=sort_AA.default, sorted_TT=sort_AT.default, sorted_TC=sort_AC.default, sorted_TG=sort_AG.default,
      sorted_CA=sort_AA.default, sorted_CT=sort_AT.default, sorted_CC=sort_AC.default, sorted_CG=sort_AG.default,
      sorted_GA=sort_AA.default, sorted_GT=sort_AT.default, sorted_GC=sort_AC.default, sorted_GG=sort_AG.default,
      sorted_null=sort_null.default,
      total_AA=sort_AA.total_reads, total_AT=sort_AT.total_reads, total_AC=sort_AC.total_reads, total_AG=sort_AG.total_reads,
      total_TA=sort_AA.total_reads, total_TT=sort_AT.total_reads, total_TC=sort_AC.total_reads, total_TG=sort_AG.total_reads,
      total_CA=sort_AA.total_reads, total_CT=sort_AT.total_reads, total_CC=sort_AC.total_reads, total_CG=sort_AG.total_reads,
      total_GA=sort_AA.total_reads, total_GT=sort_AT.total_reads, total_GC=sort_AC.total_reads, total_GG=sort_AG.total_reads,
      total_null=sort_null.total_reads
  }
  

  call mark_duplicates_split {
    input: input_bam=sort_by_bc_join.default
  }
  
  scatter(split_file in mark_duplicates_split.split_files) {
    call mark_duplicates_main {
      input: input_bam=sort_by_bc_join.default, split_file=split_file, align=extract_reads_join.out_align
    }
  }
  
  call mark_duplicates_join {
    input: output_output=mark_duplicates_main.output_bam, chunked_reporter_output=mark_duplicates_main.chunked_reporter
  }

  call count_genes_split {
    input: bam_inputs=mark_duplicates_main.output_bam, barcode_whitelist=setup_chunks.barcode_whitelist, gem_groups=extract_reads_join.gem_groups,
      reference_path=reference_path, barcode_summary=attach_bcs_and_umis_join.barcode_summary, chemistry_def=setup_chunks.chemistry_def
  }

  scatter(split_file in count_genes_split.split_files) {
    call count_genes_main {
    input: bam_inputs=mark_duplicates_main.output_bam, barcode_whitelist=setup_chunks.barcode_whitelist, gem_groups=extract_reads_join.gem_groups,
      reference_path=reference_path, barcode_summary=attach_bcs_and_umis_join.barcode_summary, chemistry_def=setup_chunks.chemistry_def,
      split_file=split_file, align=extract_reads_join.out_align
    }
  }

  call count_genes_join {
    input: matrices_h5_output=count_genes_main.matrices_h5, chunked_reporter_output=count_genes_main.chunked_reporter,
      chemistry_def=setup_chunks.chemistry_def, gem_groups=extract_reads_join.gem_groups, sample_id=sample_id
  }

  call filter_barcodes_main {
    input: sample_id=sample_id, matrices_h5=count_genes_join.matrices_h5, barcode_whitelist=setup_chunks.barcode_whitelist,
      gem_groups=extract_reads_join.gem_groups, chemistry_def=setup_chunks.chemistry_def,
      raw_fastq_summary=extract_reads_join.summary, attach_bcs_summary=attach_bcs_and_umis_join.summary,
      barcode_summary=attach_bcs_and_umis_join.barcode_summary
  }
}
