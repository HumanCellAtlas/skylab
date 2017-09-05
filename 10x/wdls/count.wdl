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
  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"

  command {

    # The sample_def maps contain a key "read_path" that refers to a path where SETUP_CHUNKS is
    # supposed to find the input fastq files. If that path is hidden within the Map, cromwell
    # will never know about it, so we have to expose it via a separate input, read_path. But,
    # the path in the Map is going to be different from the path in the Array[File] because the
    # latter is going to have a bunch of cromwell stuff prepended to it. So we need to put those
    # in a location that SETUP_CHUNKS is expecting.
    set -x
    for read_path in '${sep="\' \'" read_paths}'; do
      mkdir -p read_paths/${dollar}${lbrace}read_path##*call-setup_chunks/inputs${rbrace}
      cp -al "${dollar}read_path" read_paths/${dollar}${lbrace}read_path##*call-setup_chunks/inputs${rbrace}
    done

    /opt/tenx/run_in_10x_env.bash martian stage run SETUP_CHUNKS main \
      --sample_id '${sample_id}' \
      --chemistry_name '${chemistry_name}' \
      --sample_def '${sample_def}'
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
  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"

  command {
    set -x

    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS split \
      --chunks '${chunks}' \
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

  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"

  command {
    set -x
    
    # Move reads paths to $PWD/read_paths
    for read_path in '${sep="\' \'" read_paths}'; do
      mkdir -p read_paths/${dollar}${lbrace}read_path##*call-chunk_reads_main/shard-*/inputs${rbrace}
      cp -al "${dollar}read_path" read_paths/${dollar}${lbrace}read_path##*call-chunk_reads_main/shard-*/inputs${rbrace}
    done


    # This stage is going to call an executable that writes to a "files" directory    
    mkdir files
    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS main \
      --chunks '${chunks}' \
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

  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"

  command {
    set -x
    jq -s '.' '${sep="\' \'" out_chunks_output}' > merged_out_chunks
    
    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS join \
      --out_chunks_output merged_out_chunks
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
  command {

    /opt/tenx/run_in_10x_env.bash martian stage run EXTRACT_READS split \
      --chunks '${chunks}' \
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
 
  String dollar = '$' 
  String lbrace = '{'
  String rbrace = '}'
 
  command {
    mkdir files
    for fastq_chunk in ${sep=" " fastq_chunks}; do
      ln "${dollar}fastq_chunk" files/"${dollar}(basename ${dollar}fastq_chunk)"
    done

    /opt/tenx/run_in_10x_env.bash martian stage run EXTRACT_READS main \
      --chunks '${chunks}' \
      --chemistry_def '${chemistry_def}' \
      --barcode_whitelist ${barcode_whitelist} \
      --reads_per_file '${reads_per_file}' \
      --subsample_rate '${subsample_rate}' \
      --skip_metrics no \
      --primers '[${sep="," primers}]' \
      --split_file '${split_file}'

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
    Array[File] read2s = glob("read2s.fastq/*")
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
  Array[Array[File]]? read2s_output
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

    idx=0
    reads_arr_arg=[]
    for reads_arr in ${d}(echo '[${sep="," reads_output}]' | jq -c '.[]'); do
      mkdir -p out_reads/${d}idx
      reads_arg=[]
      for reads_path in ${d}(echo ${d}reads_arr | jq -r '.[]'); do
        new_path=out_reads/${d}idx/${d}(basename ${d}reads_path)
        ln ${d}reads_path ${d}new_path
        reads_arg=${d}( echo ${d}reads_arg | jq -c --arg v ${d}new_path '. |= .+ [$v]')
      done
      reads_arr_arg=${d}( echo ${d}reads_arr_arg | jq -c --arg v ${d}reads_arg '. |= .+ [$v | fromjson]')
      idx=${d}((idx+1))
    done
      
    idx=0
    read2s_arr_arg=[]
    for read2s_arr in ${d}(echo '[${sep="," read2s_output}]' | jq -c '.[]'); do
      mkdir -p out_read2s/${d}idx
      read2s_arg=[]
      for read2s_path in ${d}(echo ${d}read2s_arr | jq -r '.[]'); do
        new_path=out_read2s/${d}idx/${d}(basename ${d}read2s_path)
        ln ${d}read2s_path ${d}new_path
        read2s_arg=${d}( echo ${d}read2s_arg | jq -c --arg v ${d}new_path '. |= .+ [$v]')
      done
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

    gem_groups=["${d}(jq -s '.[].gem_group | tostring' ${sep=' ' split_files} | jq -rs 'join(",")')"]
    
    /opt/tenx/run_in_10x_env.bash martian stage run EXTRACT_READS join \
      --align '${align}' \
      --barcode_whitelist ${barcode_whitelist} \
      --chemistry_def '${chemistry_def}' \
      --gem_group_split "${d}gem_groups" \
      --reads_output "${d}reads_arr_arg" \
      --read2s_output "${d}read2s_arr_arg" \
      --trimmed_seqs_output "${d}trimmed_seqs_arr_arg" \
      --bam_comments_output '[${sep=',' bam_comments_output}]' \
      --gem_groups_output '[${sep=',' gem_groups_output}]' \
      --read_groups_output '[${sep=',' read_groups_output}]' \
      --chunked_reporter_output '["${sep='","' chunked_reporter_output}"]' \
      --barcode_counts_output '["${sep='","' barcode_counts_output}"]'

    jq -r '.[]' read_groups > read_groups.lines
    jq -r '.[]' gem_groups > gem_groups.lines
    jq -r '.[]' bam_comments > bam_comments.lines
    jq -r '.[]' reads > reads.lines
    jq -r '.[]' read2s > read2s.lines
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
    Array[String] read2s = read_lines("read2s.lines")
    Array[String] trimmed_seqs = read_lines("trimmed_seqs.lines")
    File out_align = "align"
  }
}

task align_reads_split {
  Array[String] read1s
  Array[String] read2s
  Array[String] read_groups
  File reference_path
    
  String d = "$"

  command {
    /opt/tenx/run_in_10x_env.bash martian stage run ALIGN_READS split \
      --reads '["${sep='","' read1s}"]' \
      --read2s '["${sep='","' read2s}"]' \
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
  File read_chunk
  File? read2_chunk
  File split_file
    
  String d = "$"

  command {

    read_group=${d}(jq -r '.read_group' ${split_file})

    /opt/tenx/run_in_10x_env.bash martian stage run ALIGN_READS main \
      --threads ${d}(nproc) \
      --read_chunk '${read_chunk}' \
      ${'-read2_chunk ' + read2_chunk} \
      --reference_path '${reference_path}' \
      --max_hits_per_read ${max_hits_per_read} \
      --reads_group "${d}read_group"
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File genome_output = "genome_output.bam"
  }
}

task attach_bcs_and_umis_split {
  String barcode_whitelist
  Array[Int] gem_groups
  Array[File] genome_inputs
  Array[File] trimmed_inputs
  
  command {
    /opt/tenx/run_in_10x_env.bash martian stage run ATTACH_BCS_AND_UMIS split \
     --gem_groups [${sep=',' gem_groups] \
     --barcode_whitelist ${barcode_whitelist} \
     --genome_inputs '["${sep='","' genome_inputs}"]' \
     --trimmed_inputs '["${sep='","' genome_inputs}"]'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] split_files = glob("martian_split_*")
  }

}

workflow count {

  String sample_id
  File sample_def
  Array[File] read_paths
  Int reads_per_file
  Float subsample_rate
  String align
  Array[Map[String, String]] primers
  File reference_path

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
    reads_output=extract_reads_main.reads, read2s_output=extract_reads_main.read2s,
    trimmed_seqs_output=extract_reads_main.trimmed_seqs, bam_comments_output=extract_reads_main.bam_comments,
    gem_groups_output=extract_reads_main.gem_groups, read_groups_output=extract_reads_main.read_groups,
    chunked_reporter_output=extract_reads_main.chunked_reporter, barcode_counts_output=extract_reads_main.barcode_counts,
    split_files=extract_reads_split.split_files
}

call align_reads_split {
  input: read1s=extract_reads_join.reads, read2s=extract_reads_join.read2s, read_groups=extract_reads_join.read_groups,
    reference_path=reference_path
}

scatter(split_pair in zip(align_reads_split, extract_reads_main.reads)) {
  call align_reads_main {
    input: reference_path=reference_path, split_file=split_pair.left, read_chunk=split_pair.right,
      max_hits_per_read=-1
  }
}
}
