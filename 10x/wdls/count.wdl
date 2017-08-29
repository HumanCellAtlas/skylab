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
  Array[File] sample_def
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

    # The sample defs are supposed to be Maps. But they're too complex for cromwell, so we input them
    # as files and read them. Python's json.load can work with them with the martian_cli.
    for sample_def_file in '${sep="\' \'" sample_def}'; do
      sample_def_strings+=("${dollar}(<${dollar}sample_def_file)")
    done

    /opt/tenx/run_in_10x_env.bash martian stage run SETUP_CHUNKS main \
      --sample_id '${sample_id}' \
      --chemistry_name '${chemistry_name}' \
      --sample_def "${dollar}${lbrace}sample_def_strings[@]${rbrace}"
  }
  

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File chunks = "chunks" 
    Map[String, String] chemistry_def = read_map("chemistry_def")
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
    echo "${dollar}(<${chunks})" | jq -c '.[]' > chunks.tmp
    # There are sometimes single quotes in strings in the chunks jsons.
    sed -i "s/'//g" chunks.tmp
    while read chunk; do
      chunks_string+=("${dollar}chunk")
    done < chunks.tmp
   

    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS split \
      --chunks "${dollar}${lbrace}chunks_string[@]${rbrace}" \
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

    for read_path in '${sep="\' \'" read_paths}'; do
      mkdir -p read_paths/${dollar}${lbrace}read_path##*call-chunk_reads_main/shard-*/inputs${rbrace}
      cp -al "${dollar}read_path" read_paths/${dollar}${lbrace}read_path##*call-chunk_reads_main/shard-*/inputs${rbrace}
    done

    chunks_string="${dollar}(<${chunks})"

    # This stage is going to call an executable that writes to a "files" directory    
    mkdir files
    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS main \
      --chunks "${dollar}chunks_string" \
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
    for out_chunks in '${sep="\' \'" out_chunks_output}'; do
      jq -c '.[]' "${dollar}out_chunks" >> out_chunks.lines
    done
    
    while read line; do
      out_chunks_array+=("[${dollar}line]")
    done < out_chunks.lines
    
    /opt/tenx/run_in_10x_env.bash martian stage run CHUNK_READS join \
      --out_chunks_output "${dollar}${lbrace}out_chunks_array[@]${rbrace}"
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
  Array[File] sample_def
  Array[File] read_paths
  Int reads_per_file

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
}
