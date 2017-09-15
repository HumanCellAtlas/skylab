# NB: This WDL use jq a lot. The manual for jq is here: https://stedolan.github.io/jq/manual/
# And you can experiment with different jq filters here: https://jqplay.org/

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/common/setup_chunks/__init__.py
task setup_chunks {
  
  # A string that will be used to name files and directories.
  String sample_id
  # See the adjacent sample_def.json.
  File sample_def
  # input files
  Array[File] r1
  Array[File] r2
  Array[File] i1
  # auto seems to be the only reasonable value
  String chemistry_name = "auto"

  command <<<

    # move all input files into read_path for cellranger
    mkdir read_path
    for fastq_file in ${sep=" " r1} ${sep=" " r2} ${sep=" " i1}; do
        mv $fastq_file read_path/
    done

    # Create the args json with the arguments
    echo "{}" > _args
    echo "$(jq --arg var "$(<${sample_def})" '. + {"sample_def": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "${chemistry_name}" '. + {"chemistry_name": $var}' _args)" > _args
    echo "$(jq --arg var "${sample_id}" '. + {"sample_id": $var}' _args)" > _args
    echo "$(jq --arg var null '. + {"custom_chemistry_def": $var | fromjson}' _args)" > _args

    # Create the outs json
    echo '{"chunks": null, "chemistry_def": null, "barcode_whitelist": null}' > _outs
    
    # Run the stage with martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/common/setup_chunks
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever
    
    # WDL-ize the outputs from the updated _outs json
    jq '.chunks' _outs > chunks.json
    jq '.chemistry_def' _outs > chemistry_def.json
    jq -r '.barcode_whitelist' _outs > barcode_whitelist.string
  >>>
  
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "2 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File chunks_setup = "chunks.json"
    File chemistry_def = "chemistry_def.json"
    String barcode_whitelist = read_string("barcode_whitelist.string")
  }
}

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/common/chunk_reads/__init__.py
task chunk_reads_split {

  # The chunks object produced by the "setup_chunks" step. This is a json that
  # defines things like read_group, gem_group, and paths to the associated read files
  File chunks_setup
  # Maximum number of reads in each chunked reads files
  Int reads_per_file

  command <<<
   	
    # Create the _args 
    echo "{}" > _args
    echo "$(jq --arg var "$(<${chunks_setup})" '. + {"chunks": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var ${reads_per_file} '. + {"reads_per_file": $var | fromjson}' _args)" > _args
    
    # Run the stage via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/common/chunk_reads
    stage_phase=split
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever
    
    # Break the chunks into separate files so we can scatter over them
    jq -c '.chunks[]' _stage_defs | split -l 1 - partition_

  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File stage_defs = "_stage_defs"
    Array[File] chunks = glob("partition_*")
  }
}

task chunk_reads_main {
  # Chunk definitions from setup_chunks
  File chunks 
  # Maximum number of reads per chunked fastq
  Int reads_per_file
  # The chunk from chunk_reads_split
  File chunk
  # input files
  Array[File] r1
  Array[File] r2
  Array[File] i1

  command <<<

    # move all input files into read_path for cellranger
    mkdir read_path
    for fastq_file in ${sep=" " r1} ${sep=" " r2} ${sep=" " i1}; do
        mv $fastq_file read_path/
    done

    # This stage is going to call an executable that writes to a "files" directory
    mkdir files

    # Create the args object
    echo "{}" > _args
    echo "$(jq --arg var "$(<${chunks})" '. + {"chunks": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var ${reads_per_file} '. + {"reads_per_file": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "$(<${chunk})" '. |=  .+ ($var | fromjson) ' _args)" > _args
   
    # Create the outs object
    echo '{"out_chunks": null}' > _outs
    
    # Execute via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/common/chunk_reads
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . files run_file_whatever
    
    # Read the out_chunks definitions into a file
    jq '.out_chunks' _outs > out_chunks.json

  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File outs = "_outs"
    Array[File] fastq_chunks = glob("files/*.fastq")
  }
}

task chunk_reads_join {
  # Array of outs objects from the chunk_reads_main steps. All this step is doing
  # is joining some dicts together.
  Array[File] outs

  command <<<
    
    # Create the args
    echo "{}" > _args
    ehco "{}" > _chunk_defs
    jq -s '.' '${sep="\' \'" outs}' > _chunk_outs
    echo '{"out_chunks": null}' > _outs
    
    # Execute via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/common/chunk_reads
    stage_phase=join
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever
    
    # Write out_chunks to their own file
    jq '.out_chunks' _outs > out_chunks.json
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File out_chunks = "out_chunks.json"
  }
}

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/common/extract_reads/__init__.py
task extract_reads_split {
  # The chunks defined by chunk_reads
  File out_chunks
  # String for the barcode whitelist, usually "737K-august-2016"
  String barcode_whitelist
  
  command <<<
    # Create the _args 
    echo '{"initial_reads": null}' > _args
    echo "$(jq --arg var "$(<${out_chunks})" '. + {"chunks": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "${barcode_whitelist}" '. + {"barcode_whitelist": $var}' _args)" > _args
    
    # Run via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/common/extract_reads
    stage_phase=split
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever

    # Break the chunks into separate files so we can scatter over them
    jq -c '.chunks[]' _stage_defs | split -l 1 - partition_
  >>>
 
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File stage_defs = "_stage_defs"
    Array[File] chunks = glob("partition_*")
  }
}

task extract_reads_main {
  # out_chunks definitions from chunk_reads
  File out_chunks
  # chemistry json produced by setup_chunks
  File chemistry_def
  # String for the barcode whitelist, usually "737K-august-2016"
  String barcode_whitelist
  # Maximum number of reads per file
  Int reads_per_file
  Float subsample_rate
  # Primers defined in the input to the workflow
  Array[Map[String, String]] primers

  Array[File] fastq_chunks
  File chunk
 
  command <<<
    
    # Move the fastq chuks to cwd/files. That's how they're referenced in the out_chunks json
    mkdir files
    for fastq_chunk in ${sep=" " fastq_chunks}; do
      ln "$fastq_chunk" files/"$(basename $fastq_chunk)"
    done

    # Create the _args 
    echo '{"rna_read_length": null, "initial_read": null, "skip_metrics": false}' > _args
    echo "$(jq --arg var "$(<${out_chunks})" '. + {"chunks": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "$(<${chemistry_def})" '. + {"chemistry_def": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "${barcode_whitelist}" '. + {"barcode_whitelist": $var}' _args)" > _args
    echo "$(jq --arg var ${reads_per_file} '. + {"reads_per_file": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var ${subsample_rate} '. + {"subsample_rate": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var '[${sep="," primers}]' '. + {"primers": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "$(<${chunk})" '. |=  .+ ($var | fromjson) ' _args)" > _args
    
    # Create the _outs
    echo '{"chunked_reporter": "chunked_reporter.pickle", "summary": "summary.json", "barcode_counts": "barcode_counts.json", ' > _outs
    echo '"reads": "reads", "read2s": "read2s", "trimmed_seqs": "trimmed_seqs", "gem_groups": null, ' >> _outs
    echo '"read_groups": null, "align": null, "bam_comments": null}' >> _outs
    
    # Run via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/common/extract_reads
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . files run_file_whatever
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    Array[File] trimmed_seqs = glob("trimmed_seqs/*")
    Array[File] reads = glob("reads/*")
    File barcode_counts = "barcode_counts.json"
    File chunked_reporter = "chunked_reporter.pickle"
    File outs = "_outs"
  }
}

task extract_reads_join {
  # String that defines parameters to the aligner. An input to the workflow
  String align 
  String barcode_whitelist
  File chemistry_def
  # Output jsons produced by the extract_reads_main steps
  Array[File] outs
  # Output of extract_reads_split
  File stage_defs

  # Files from the main steps. We have to be explicit about these because
  # Martian doesn't handle file staging
  Array[Array[File]] reads_from_mains
  Array[Array[File]] trimmed_seqs_from_mains
  Array[File] chunked_reporter_from_mains
  Array[File] barcode_counts_from_mains
    
  # These are defined so we can work with bash arrays with it being interpreted
  # as WDL
  String l = "{"
  String r = "}"

  command <<<

    # Create args
    echo '{}' > _args
    echo "$(jq --arg var "${barcode_whitelist}" '. + {"barcode_whitelist": $var}' _args)" > _args
    echo "$(jq --arg var '${align}' '. + {"align": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "$(<${chemistry_def})" '. + {"chemistry_def": $var | fromjson}' _args)" > _args
   
    # Create chunk_defs
    jq '.chunks' '${stage_defs}' > _chunk_defs
    
    # Create chunk_outs. The outs files produced by the main stages contain paths that point
    # to files in those stages. We need to replace those with paths valid in this task.
    reads_array=('${sep="\' \'" reads_from_mains}')
    trimmed_seqs_array=('${sep="\' \'" trimmed_seqs_from_mains}')
    chunked_reporter_array=(${sep=" " chunked_reporter_from_mains})
    barcode_counts_array=(${sep=" " barcode_counts_from_mains})
    idx=0
    for out in '${sep="\' \'" outs}'; do
      reads_path="$${l}reads_array[$idx]${r}"
      trimmed_seqs_path="$${l}trimmed_seqs_array[$idx]${r}"
      chunked_reporter_path="$${l}chunked_reporter_array[$idx]${r}"
      barcode_counts_path="$${l}barcode_counts_array[$idx]${r}"
      echo "$(jq --arg var "$reads_path" '.reads = ( $var | fromjson )' $out)" > $out
      echo "$(jq --arg var "$trimmed_seqs_path" '.trimmed_seqs = ( $var | fromjson )' $out)" > $out
      echo "$(jq --arg var "$chunked_reporter_path" '.chunked_reporter = $var' $out)" > $out
      echo "$(jq --arg var "$barcode_counts_path" '.barcode_counts = $var' $out)" > $out
      idx=$((idx+1))
    done
    jq -s '.' '${sep="\' \'" outs}' > _chunk_outs
  
    # Create outs
    echo '{"summary": "summary.json", "barcode_counts": "barcode_counts.json", "reads": null, ' > _outs
    echo '"reads2": null, "trimmed_seqs": "trimmed_seqs", "gem_groups": null, "read_groups": null, ' >> _outs
    echo '"align": null, "bam_comments": null}' >> _outs
    
    # Run the martian stage
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/common/extract_reads
    stage_phase=join
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . files run_file_whatever

    # WDL-ize outputs
    jq -r '.read_groups[]' _outs > read_groups.lines
    jq -r '.gem_groups[]' _outs > gem_groups.lines
    jq -r '.bam_comments[]' _outs > bam_comments.lines
    
    # Output the list of reads files and trimmed_seqs files
    # First move the reads and trimmed seqs to a path in cwd so we won't
    # have to do this in future stages
    mkdir reads
    mkdir trimmed_seqs
    idx=0
    for reads_file in $(jq -r '.reads[]' _outs); do
      ln "$reads_file" reads/$idx.fastq
      idx=$((idx+1))
    done
    idx=0
    for trimmed_seqs_file in $(jq -r '.trimmed_seqs[]' _outs); do
      ln "$trimmed_seqs_file" trimmed_seqs/$idx.bam
      idx=$((idx+1))
    done
    
    # Then update the _outs with those paths
    echo "$(jq '.reads = []' _outs)" > _outs
    echo "$(jq '.trimmed_seqs = []' _outs)" > _outs

    for reads_file in reads/*; do
      echo "$(jq --arg var $reads_file '.reads |= .+ [$var]' _outs)" > _outs
    done

    for trimmed_seqs_file in trimmed_seqs/*; do
      echo "$(jq --arg var $trimmed_seqs_file '.trimmed_seqs |= .+ [$var]' _outs)" > _outs
    done
    
    # And finally write those paths to something WDL can read
    jq -r '.reads[]' _outs > read_paths.lines
    jq -r '.trimmed_seqs[]' _outs > trimmed_seq_paths.lines
    jq '.align' _outs > align.json
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File summary = "summary.json"
    File barcode_counts = "barcode_counts.json"
    File out_align = "align.json"
    Array[File] reads = glob("reads/*")
    Array[File] trimmed_seqs = glob("trimmed_seqs/*")
    Array[String] read_paths = read_lines("read_paths.lines")
    Array[String] trimmed_seq_paths = read_lines("trimmed_seq_paths.lines")
    Array[String] read_groups = read_lines("read_groups.lines")
    Array[Int] gem_groups = read_lines("gem_groups.lines")
    Array[String] bam_comments = read_lines("bam_comments.lines")
  }
}

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/align_reads/__init__.py
task align_reads_split {
  # The split phase doesn't really need the reads files, so it's fine just giving it strings
  # that are (invalid) paths. It's just splitting the list up.
  Array[String] read_paths
  Array[String] read_groups
  File reference_path
  
  command <<<
    # todo is this necessary here? is the reference itself used or just the path to it?
    mkdir genome_reference && tar -xf ${reference_path} -C genome_reference --strip-components 1

    # Create the args
    echo '{"threads": 2}' > _args
    read_paths_json='["${sep='","' read_paths}"]'
    echo "$(jq --arg var "$read_paths_json" '. + {"reads": ( $var | fromjson )}' _args)" > _args
    read_groups_json='["${sep='","' read_groups}"]'
    echo "$(jq --arg var "$read_groups_json" '. + {"read_groups": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "[]" '. + {"read2s": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "genome_reference" '. + {"reference_path": $var}' _args)" > _args
    
    # Run martian
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/align_reads
    stage_phase=split
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever

    # Break the chunks into separate files so we can scatter over them
    jq -c '.chunks[]' _stage_defs | split -l 1 - partition_
  >>>
 
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File stage_defs = "_stage_defs"
    Array[File] chunks = glob("partition_*")
  }
}

task align_reads_main {

  Int max_hits_per_read
  File reference_path
  # A fastq produced by extract_reads 
  File chunked_fastq
  # The chunk definition produced by align_reads_split
  File chunk

  command <<<

    mkdir genome_reference && tar -xf ${reference_path} -C genome_reference --strip-components 1

    # Create args
    threads=$(nproc)
    echo '{"threads": '$threads', "read2_chunk": null}'  > _args
    echo "$(jq --arg var genome_reference '. + {"reference_path": $var}' _args)" > _args
    echo "$(jq --arg var ${max_hits_per_read} '. + {"max_hits_per_read": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "$(<${chunk})" '. |=  .+ ($var | fromjson) ' _args)" > _args
    echo "$(jq --arg var "${chunked_fastq}" '. + {"read_chunk": $var}' _args)" > _args
    
    # Create outs
    echo '{"genome_output": "genome_output.bam"}' > _outs

    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/align_reads
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File genome_output = "genome_output.bam"
  }
}

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/attach_bcs_and_umis/__init__.py
task attach_bcs_and_umis_main {
  File reference_path
  String barcode_whitelist
  # json of barcode counts produced by extract_reads_join
  File barcode_counts
  File chemistry_def
  Array[Int] gem_groups
  Int umi_min_qual_threshold
  Array[String] bam_comments
  
  # See the workflow definition, but here unpack the pair structure that these
  # main stages are being scattered over
  Pair[Int, Pair[File, File]] chunk
  Int gem_group = chunk.left
  Pair[File, File] chunk_right = chunk.right
  File trimmed_seqs = chunk_right.left
  File genome_output= chunk_right.right

  command <<<

    # todo check if reference or just path is needed here.
    mkdir genome_reference && tar -xf ${reference_path} -C genome_reference --strip-components 1

    # Create args. Lots of args
    echo '{"rescue_multimappers": true, "correct_barcodes": true, "paired_end": false, "skip_metrics": false, ' > _args
    echo '"barcode_confidence_threshold": 0.975, "annotation_params": null}' >> _args
    echo "$(jq --arg var "${barcode_whitelist}" '. + {"barcode_whitelist": $var}' _args)" > _args
    echo "$(jq --arg var genome_reference '. + {"reference_path": $var}' _args)" > _args
    echo "$(jq --arg var ${barcode_counts} '. + {"barcode_counts": $var}' _args)" > _args
    echo "$(jq --arg var "$(<${chemistry_def})" '. + {"chemistry_def": $var | fromjson}' _args)" > _args
    bam_comments_json='["${sep='","' bam_comments}"]'
    echo "$(jq --arg var "$bam_comments_json" '. + {"bam_comments": ( $var | fromjson )}' _args)" > _args
    gem_groups_json='[${sep=',' gem_groups}]'
    echo "$(jq --arg var "$gem_groups_json" '. + {"gem_groups": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "${gem_group}" '. + {"gem_group": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "${umi_min_qual_threshold}" '. + {"umi_min_qual_threshold": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "${genome_output}" '. + {"chunk_genome_input": $var}' _args)" > _args
    echo "$(jq --arg var "${trimmed_seqs}" '. + {"chunk_trimmed_input": $var}' _args)" > _args
    
    # Create the outs
    echo '{"chunked_reporter": "chunked_reporter.pickle", "output": "output.bam", "num_alignments": null}' > _outs
    
    # Run via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/attach_bcs_and_umis
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever
    
    # Write num_alignments to something WDL can read
    jq -r '.num_alignments' _outs > num_alignments.int
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    Int num_alignments = read_int("num_alignments.int")
    File output_bam = "output.bam"
    File chunked_reporter = "chunked_reporter.pickle"
    File outs = "_outs"
  }
}


task attach_bcs_and_umis_join {
  File reference_path

  # Outputs from the main phases
  Array[Int] num_alignments_from_mains
  Array[File] output_from_mains
  Array[File] chunked_reporter_from_mains
  # outs jsons from the main phases
  Array[File] outs
  
  # Needed to work with bash arrays
  String l = "{"
  String r = "}"

  command <<<

    # todo check if reference or just path is needed here.
    mkdir genome_reference && tar -xf ${reference_path} -C genome_reference --strip-components 1

    # Create the args
    echo '{}' > _args
    echo "$(jq --arg var genome_reference '. + {"reference_path": $var}' _args)" > _args
  
    # Create chunk_outs. Replace paths from the main stages with paths that are valid in
    # this join stage
    num_alignments_array=(${sep=" " num_alignments_from_mains})
    output_array=('${sep="\' \'" output_from_mains}')
    chunked_reporter_array=('${sep="\' \'" chunked_reporter_from_mains}')
    idx=0
    for out in '${sep="\' \'" outs}'; do
      num_alignments=$${l}num_alignments_array[$idx]${r}
      chunked_reporter_path="$${l}chunked_reporter_array[$idx]${r}"
      output_path="$${l}output_array[$idx]${r}"
      echo "$(jq --arg var "$num_alignments" '.reads = ( $var | fromjson )' $out)" > $out
      echo "$(jq --arg var "$chunked_reporter_path" '.chunked_reporter = $var' $out)" > $out
      echo "$(jq --arg var "$output_path" '.output = $var' $out)" > $out
      idx=$((idx+1))
    done
    jq -s '.' '${sep="\' \'" outs}' > _chunk_outs

    # Create the outs
    echo '{"summary": "summary.json", "barcode_summary": "barcode_summary.h5", "num_alignments": null}' > _outs
    
    # Create a dummy _chunk_defs. Stage doesn't use it, but the file has to be there.
    echo '{}' > _chunk_defs

    # Run via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/attach_bcs_and_umis
    stage_phase=join
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever
    
    # Make num_alignments readable by WDL
    jq '.num_alignments[]' _outs > num_alignments.lines
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File summary = "summary.json"
    File barcode_summary = "barcode_summary.h5"
    Array[Int] num_alignments = read_lines("num_alignments.lines")
  }
}

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/bucket_reads_by_bc/__init__.py
task bucket_by_bc_split {
  # The BAMs produced by all the attach_umis_and_bcs main steps
  Array[File] inputs
  # The alignment counts produced by all the attach_umis_and_bcs main steps
  Array[Int] num_alignments

  command <<<
    
    # Create args
    echo '{"nbases": 2}' >  _args
    num_alignments_json='[${sep=", " num_alignments}]'
    echo "$(jq --arg var "$num_alignments_json" '. + {"num_alignments": ( $var | fromjson )}' _args)" > _args
    inputs_json='["${sep='\", \"' inputs}"]'
    echo "$(jq --arg var "$inputs_json" '. + {"inputs": ( $var | fromjson )}' _args)" > _args
    
    # Run martian
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/bucket_reads_by_bc
    stage_phase=split
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever

    # Break the chunks into separate files so we can scatter over them
    jq '.chunks[0].read_groups' _stage_defs > read_groups
  >>>
 
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File stage_defs = "_stage_defs"
    File read_groups = "read_groups"
  }
}

task bucket_by_bc_main {
  # One of the BAMs produced by attach_umis
  File chunk_input
  # json with a list of the read groups
  File read_groups
  

  command <<<
    
    # Create the args. cellranger hard codes nbases to 2
    echo '{"nbases": 2}' >  _args
    echo "$(jq --arg var "${chunk_input}" '. + {"chunk_input": $var}' _args)" > _args
    echo "$(jq --arg var "$(<${read_groups})" '. + {"read_groups": ( $var | fromjson )}' _args)" > _args

    mkdir files

    # Create the outs
    echo '{"chunked_reporter": "chunked_reporter.pickle", "output": "output.bam", "num_alignments": null}' > _outs
    
    # Run via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/bucket_reads_by_bc
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . files run_file_whatever
    
    # Pull the buckets definitions into a separate json file
    jq '.buckets' _outs > buckets
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File outs = "_outs"
    File buckets = "buckets"
    # This splits up the BAM into 17 buckets that we have to parallelize over. I'm not sure
    # how to do this concisely, so just being very verbose about it.
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

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/sort_reads_by_bc/__init__.py
task sort_by_bc_main {
  # A set of BAM from the bucket_by_bc steps. So this would be a list of "bc_AT.bam" files,
  # for example
  Array[File] bucketed_bams
  
  command <<<
    # All the files have the same name, so move the bucketed bams to cwd with unique names
    idx=0
    bucket='[]'
    for bam_path in '${sep="\' \'" bucketed_bams}'; do
      mkdir -p bucketed_bams/$idx
      new_path=bucketed_bams/$idx/"$(basename $bam_path)"
      ln "$bam_path" $new_path
      bucket=$(echo $bucket | jq -c --arg v $new_path '. |= .+ [$v]')
      idx=$((idx+1))
    done
    
    echo '{}' > _args
    echo "$(jq --arg var "$bucket" '. + {"bucket": ( $var | fromjson )}' _args)" > _args
    
    echo '{"default": "default.bam", "total_reads": null}' > _outs
    mkdir files

    # Run via martian_shell
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/sort_reads_by_bc
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . files run_file_whatever
    
    jq '.total_reads' _outs > total_reads.int
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File default = "default.bam"
    Int total_reads = read_int("total_reads.int")
    File outs = "_outs"
  }
}

task sort_by_bc_join {
  # Merge the 17 buckets together.
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

  command <<<
    set -x
    
    echo '{}' > _args
    echo '[{"prefix": "AA"}, ' >> _chunk_defs
    echo '{"prefix": "AT"}, ' >> _chunk_defs
    echo '{"prefix": "AC"}, ' >> _chunk_defs
    echo '{"prefix": "AG"}, ' >> _chunk_defs
    echo '{"prefix": "TA"}, ' >> _chunk_defs
    echo '{"prefix": "TT"}, ' >> _chunk_defs
    echo '{"prefix": "TC"}, ' >> _chunk_defs
    echo '{"prefix": "TG"}, ' >> _chunk_defs
    echo '{"prefix": "CA"}, ' >> _chunk_defs
    echo '{"prefix": "CT"}, ' >> _chunk_defs
    echo '{"prefix": "CC"}, ' >> _chunk_defs
    echo '{"prefix": "CG"}, ' >> _chunk_defs
    echo '{"prefix": "GA"}, ' >> _chunk_defs
    echo '{"prefix": "GT"}, ' >> _chunk_defs
    echo '{"prefix": "GC"}, ' >> _chunk_defs
    echo '{"prefix": "GG"}, ' >> _chunk_defs
    echo '{"prefix": ""}]' >> _chunk_defs
    
    echo '[{"default": "${sorted_AA}", "total_reads": ${total_AA}}, ' >> _chunk_outs
    echo '{"default": "${sorted_AT}", "total_reads": ${total_AT}}, ' >> _chunk_outs
    echo '{"default": "${sorted_AC}", "total_reads": ${total_AC}}, ' >> _chunk_outs
    echo '{"default": "${sorted_AG}", "total_reads": ${total_AG}}, ' >> _chunk_outs
    echo '{"default": "${sorted_TA}", "total_reads": ${total_TA}}, ' >> _chunk_outs
    echo '{"default": "${sorted_TT}", "total_reads": ${total_TT}}, ' >> _chunk_outs
    echo '{"default": "${sorted_TC}", "total_reads": ${total_TC}}, ' >> _chunk_outs
    echo '{"default": "${sorted_TG}", "total_reads": ${total_TG}}, ' >> _chunk_outs
    echo '{"default": "${sorted_CA}", "total_reads": ${total_CA}}, ' >> _chunk_outs
    echo '{"default": "${sorted_CT}", "total_reads": ${total_CT}}, ' >> _chunk_outs
    echo '{"default": "${sorted_CC}", "total_reads": ${total_CC}}, ' >> _chunk_outs
    echo '{"default": "${sorted_CG}", "total_reads": ${total_CG}}, ' >> _chunk_outs
    echo '{"default": "${sorted_GA}", "total_reads": ${total_GA}}, ' >> _chunk_outs
    echo '{"default": "${sorted_GT}", "total_reads": ${total_GT}}, ' >> _chunk_outs
    echo '{"default": "${sorted_GC}", "total_reads": ${total_GC}}, ' >> _chunk_outs
    echo '{"default": "${sorted_GG}", "total_reads": ${total_GG}}, ' >> _chunk_outs
    echo '{"default": "${sorted_null}", "total_reads": ${total_null}}]' >> _chunk_outs
    
    echo '{"default": "default.bam", "total_reads": null}' > _outs

    mkdir files
    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/sort_reads_by_bc
    stage_phase=join
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . files run_file_whatever
    
    jq '.total_reads' _outs > total_reads.int
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File default = "default.bam"
    Int total_reads = read_int("total_reads.int")
  }
}

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/mark_duplicates/__init__.py
task mark_duplicates_split {
  File input_bam
  
  command <<<
    # Create args
    echo '{"mem_gb": null}' >  _args
    echo "$(jq --arg var "${input_bam}" '. + {"input": $var}' _args)" > _args

    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/mark_duplicates
    stage_phase=split
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever

    # Break the chunks into separate files so we can scatter over them
    jq -c '.chunks[]' _stage_defs | split -l 1 - partition_
  >>>
 
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File stage_defs = "_stage_defs"
    Array[File] chunks = glob("partition_*")
  }
}

task mark_duplicates_main {
  File input_bam
  File align
  File reference_path
  File chunk

  command <<<

    # todo check if reference or just path is needed here.
    mkdir genome_reference && tar -xf ${reference_path} -C genome_reference --strip-components 1

    echo '{}'  > _args
    echo "$(jq --arg var ${input_bam} '. + {"input": $var}' _args)" > _args
    echo "$(jq --arg var genome_reference '. + {"reference_path": $var}' _args)" > _args
    echo "$(jq --arg var "$(<${align})" '. + {"align": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "$(<${chunk})" '. |=  .+ ($var | fromjson) ' _args)" > _args
    
    echo '{"chunked_reporter": "chunked_reporter.pickle", "output": "output.bam", "summary": "summary.json"}' > _outs

    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/mark_duplicates
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever
  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File chunked_reporter = "chunked_reporter.pickle"
    File output_bam = "output.bam"
    File out = "_outs"
  }
}

task mark_duplicates_join {
  Array[File] output_bam_from_mains
  Array[File] chunked_reporter_from_mains
  Array[File] outs
  
  String l = "{"
  String r = "}"
  
  command <<<

    echo '{}' > _args
    echo '{}' > _chunk_defs

    echo '{"chunked_reporter": "chunked_reporter.pickle", "output": "output.bam", "summary": "summary.json"}' > _outs

    # Create chunk_outs
    bam_array=('${sep="\' \'" output_bam_from_mains}')
    chunked_reporter_array=(${sep=" " chunked_reporter_from_mains})
    idx=0
    for out in '${sep="\' \'" outs}'; do
      bam_path="$${l}reads_array[$idx]${r}"
      chunked_reporter_path="$${l}chunked_reporter_array[$idx]${r}"
      echo "$(jq --arg var "$bam_path" '.bam = $var ' $out)" > $out
      echo "$(jq --arg var "$chunked_reporter_path" '.chunked_reporter = $var' $out)" > $out
      idx=$((idx+1))
    done
    jq -s '.' '${sep="\' \'" outs}' > _chunk_outs

    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/mark_duplicates
    stage_phase=join
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever
  >>>
  
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File summary = "summary.json"
  }
}

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/count_genes/__init__.py
task count_genes_split {
  Array[File] bam_inputs
  String barcode_whitelist
  Array[Int] gem_groups
  File reference_path
  File barcode_summary
  File chemistry_def


  command <<<

    # todo check if reference or just path is needed here.
    mkdir genome_reference && tar -xf ${reference_path} -C genome_reference --strip-components 1

    idx=0
    linked_bam_inputs='[]'
    for bam_path in '${sep="\' \'" bam_inputs}'; do
      mkdir -p bam_inputs/$idx
      new_path=bam_inputs/$idx/"$(basename $bam_path)"
      ln "$bam_path" $new_path
      linked_bam_inputs=$(echo $linked_bam_inputs | jq -c --arg v $new_path '. |= .+ [$v]')
      idx=$((idx+1))
    done

    # Create args
    echo '{}' > _args
    echo "$(jq --arg var "${barcode_whitelist}" '. + {"barcode_whitelist": $var}' _args)" > _args
    gem_groups_json='[${sep=',' gem_groups}]'
    echo "$(jq --arg var "$gem_groups_json" '. + {"gem_groups": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "genome_reference" '. + {"reference_path": $var}' _args)" > _args
    echo "$(jq --arg var "${barcode_summary}" '. + {"barcode_summary": $var}' _args)" > _args
    echo "$(jq --arg var "$(<${chemistry_def})" '. + {"chemistry_def": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "$linked_bam_inputs" '. + {"inputs": ($var | fromjson)}' _args)" > _args

    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/count_genes
    stage_phase=split
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever

    # Break the chunks into separate files so we can scatter over them
    jq -c '.chunks[]' _stage_defs | split -l 1 - partition_
  >>>
 
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File stage_defs = "_stage_defs"
    Array[File] chunks = glob("partition_*")
  }
}


task count_genes_main {
  Array[File] bam_inputs
  String barcode_whitelist
  Array[Int] gem_groups
  File reference_path
  File barcode_summary
  File chemistry_def
  File chunk
  File align

  command <<<

    # todo check if reference or just path is needed here.
    mkdir genome_reference && tar -xf ${reference_path} -C genome_reference --strip-components 1

    idx=0
    linked_bam_inputs='[]'
    for bam_path in '${sep="\' \'" bam_inputs}'; do
      mkdir -p bam_inputs/$idx
      new_path=bam_inputs/$idx/"$(basename $bam_path)"
      ln "$bam_path" $new_path
      linked_bam_inputs=$(echo $linked_bam_inputs | jq -c --arg v $new_path '. |= .+ [$v]')
      idx=$((idx+1))
    done

    # Create args
    echo '{}' > _args
    echo "$(jq --arg var "${barcode_whitelist}" '. + {"barcode_whitelist": $var}' _args)" > _args
    gem_groups_json='[${sep=',' gem_groups}]'
    echo "$(jq --arg var "$gem_groups_json" '. + {"gem_groups": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "genome_reference" '. + {"reference_path": $var}' _args)" > _args
    echo "$(jq --arg var "${barcode_summary}" '. + {"barcode_summary": $var}' _args)" > _args
    echo "$(jq --arg var "$(<${chemistry_def})" '. + {"chemistry_def": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "$(<${align})" '. + {"align": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "$linked_bam_inputs" '. + {"inputs": ($var | fromjson)}' _args)" > _args
    echo "$(jq --arg var "$(<${chunk})" '. |=  .+ ($var | fromjson) ' _args)" > _args
    
    echo '{"chunked_reporter": "chunked_reporter.pickle", "matrices_h5": "matrices.h5"}' > _outs

    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/count_genes
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever

  >>>

  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File matrices_h5 = "matrices.h5"
    File chunked_reporter = "chunked_reporter.pickle"
    File outs = "_outs"
  }
}

task count_genes_join {
  Array[Int] gem_groups
  File chemistry_def
  String sample_id

  Array[File] outs
  Array[File] matrices_h5_from_mains
  Array[File] chunked_reporter_from_mains
  
  String l = "{"
  String r = "}"

  command <<<

    # Create args
    echo '{}' > _args
    gem_groups_json='[${sep=',' gem_groups}]'
    echo "$(jq --arg var "$(<${chemistry_def})" '. + {"chemistry_def": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "$gem_groups_json" '. + {"gem_groups": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "${sample_id}" '. + {"sample_id": $var}' _args)" > _args
    
    echo '{}' > _chunk_defs

    # Create the chunk outs
    matrices_array=(${sep=" " matrices_h5_from_mains})
    chunked_reporter_array=(${sep=" " chunked_reporter_from_mains})
    idx=0
    for out in '${sep="\' \'" outs}'; do
      matrices_path="$${l}matrices_array[$idx]${r}"
      chunked_reporter_path="$${l}chunked_reporter_array[$idx]${r}"
      echo "$(jq --arg var "$matrices_path" '.matrices_h5 = $var ' $out)" > $out
      echo "$(jq --arg var "$chunked_reporter_path" '.chunked_reporter = $var' $out)" > $out
      idx=$((idx+1))
    done
    jq -s '.' '${sep="\' \'" outs}' > _chunk_outs
    
    echo '{"barcode_summary": "barcode_summary.h5", "matrices_h5": "matrices.h5", "reporter_summary": "reporter_summary.json", "matrices_mex": "matrices_mex"}' > _outs

    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/count_genes
    stage_phase=join
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever

    # tar the matrices_mex folder so it can be extracted from the task
    tar -czf matrices_mex.tar.gz matrices_mex
  >>>
  
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File matrices_h5 = "matrices.h5"
    File matrices_mex = "matrices_mex.tar.gz"
    File reporter_summary = "reporter_summary.json"
    File barcode_summary = "barcode_summary.h5"
  }

}

# See https://github.com/10XGenomics/cellranger/blob/master/mro/stages/counter/filter_barcodes/__init__.py
task filter_barcodes {

  String sample_id
  File matrices_h5
  File raw_fastq_summary
  File attach_bcs_summary
  String barcode_whitelist
  Array[Int] gem_groups
  File chemistry_def
  File barcode_summary


  command <<<
    # Create args
    echo '{"force_cells": null, "recovered_cells": null, "cell_barcodes": null}' > _args
    echo "$(jq --arg var "${sample_id}" '. + {"sample_id": $var}' _args)" > _args
    echo "$(jq --arg var "${matrices_h5}" '. + {"matrices_h5": $var}' _args)" > _args
    echo "$(jq --arg var "${raw_fastq_summary}" '. + {"raw_fastq_summary": $var}' _args)" > _args
    echo "$(jq --arg var "${attach_bcs_summary}" '. + {"attach_bcs_summary": $var}' _args)" > _args
    echo "$(jq --arg var "${barcode_whitelist}" '. + {"barcode_whitelist": $var}' _args)" > _args
    gem_groups_json='[${sep=',' gem_groups}]'
    echo "$(jq --arg var "$gem_groups_json" '. + {"gem_groups": ( $var | fromjson )}' _args)" > _args
    echo "$(jq --arg var "$(<${chemistry_def})" '. + {"chemistry_def": $var | fromjson}' _args)" > _args
    echo "$(jq --arg var "${barcode_summary}" '. + {"barcode_summary": $var}' _args)" > _args
    

    # Create outs
    echo '{"summary": "summary.json", "filtered_barcodes": "filtered_barcodes.csv", ' >> _outs
    echo '"filtered_matrices_h5": "filtered_matrices.h5", "filtered_matrices_mex": "filtered_matrices_mex"}' >> _outs

    mv /_jobinfo .  
    stage_path=/cellranger/mro/stages/counter/filter_barcodes
    stage_phase=main
    /martian/adapters/python/martian_shell.py $stage_path $stage_phase . . run_file_whatever

    # tar the filtered_matrices_mex folder so it can be extracted from the task
    tar -czf filtered_matrices_mex.tar.gz filtered_matrices_mex

  >>>
  
  runtime {
    docker: "marcusczi/cellranger_clean:cromwell"
    memory: "30 GB"
    disks: "local-disk 50 HDD"
  }

  output {
    File summary = "summary.json"
    File filtered_barcodes = "filtered_barcodes.csv"
    File filtered_matrices_h5 = "filtered_matrices.h5"
    File filtered_matrices_mex = "filtered_matrices_mex.tar.gz"
  }
}

workflow count {

  File sample_def
  Array[File] r1
  Array[File] r2
  Array[File] i1
  String sample_id
  Int reads_per_file
  Float subsample_rate
  Array[Map[String, String]] primers
  String align
  File reference_path
  Int umi_min_qual_threshold

  call setup_chunks {
    input:
      sample_id = sample_id,
      sample_def = sample_def,
      r1 = r1,
      r2 = r2,
      i1 = i1
  }

  call chunk_reads_split {
    input:
      chunks_setup = setup_chunks.chunks_setup,
      reads_per_file = reads_per_file
  }

  scatter(chunk in chunk_reads_split.chunks) {
    call chunk_reads_main {
      input:
        chunk = chunk,
        chunks = setup_chunks.chunks_setup,
        reads_per_file = reads_per_file,
        r1 = r1,
        r2 = r2,
        i1 = i1
    }
  }

  call chunk_reads_join {
    input:
      outs = chunk_reads_main.outs
  }

  call extract_reads_split {
    input:
      out_chunks = chunk_reads_join.out_chunks,
      barcode_whitelist = setup_chunks.barcode_whitelist
  }

  scatter(split_pair in zip(extract_reads_split.chunks, chunk_reads_main.fastq_chunks)) {
      call extract_reads_main {
        input:
          out_chunks = chunk_reads_join.out_chunks,
          barcode_whitelist = setup_chunks.barcode_whitelist,
          primers = primers,
          reads_per_file = reads_per_file,
          subsample_rate = subsample_rate,
          chemistry_def = setup_chunks.chemistry_def,
          chunk = split_pair.left,
          fastq_chunks = split_pair.right
    }
  }

  call extract_reads_join {
    input:
      align = align,
      barcode_whitelist = setup_chunks.barcode_whitelist,
      chemistry_def = setup_chunks.chemistry_def,
      outs = extract_reads_main.outs,
      stage_defs = extract_reads_split.stage_defs,
      reads_from_mains = extract_reads_main.reads,
      trimmed_seqs_from_mains = extract_reads_main.trimmed_seqs,
      chunked_reporter_from_mains = extract_reads_main.chunked_reporter,
      barcode_counts_from_mains = extract_reads_main.barcode_counts
  }

  call align_reads_split {
    input:
      read_groups = extract_reads_join.read_groups,
      read_paths = extract_reads_join.read_paths,
      reference_path = reference_path
  }

  scatter(split_pair in zip(align_reads_split.chunks, extract_reads_join.reads)) {
    call align_reads_main {
      input:
        reference_path = reference_path,
        chunk = split_pair.left,
        chunked_fastq = split_pair.right,
        max_hits_per_read = -1
    }
  }

  Array[Pair[File, File]] trimmed_and_aligned  =  zip(extract_reads_join.trimmed_seqs, align_reads_main.genome_output)
  Array[Pair[Int, Pair[File, File]]] gem_group_trimmed_and_aligned  =  zip(extract_reads_join.gem_groups, trimmed_and_aligned)

  scatter(chunk in gem_group_trimmed_and_aligned) {
    call attach_bcs_and_umis_main {
      input:
        reference_path = reference_path,
        barcode_whitelist = setup_chunks.barcode_whitelist,
        chemistry_def = setup_chunks.chemistry_def,
        gem_groups = extract_reads_join.gem_groups,
        umi_min_qual_threshold = umi_min_qual_threshold,
        barcode_counts = extract_reads_join.barcode_counts,
        chunk = chunk,
        bam_comments = extract_reads_join.bam_comments 
    }
  }

  call attach_bcs_and_umis_join {
    input:
      reference_path = reference_path,
      num_alignments_from_mains = attach_bcs_and_umis_main.num_alignments,
      output_from_mains = attach_bcs_and_umis_main.output_bam,
      chunked_reporter_from_mains = attach_bcs_and_umis_main.chunked_reporter,
      outs = attach_bcs_and_umis_main.outs
  }

  call bucket_by_bc_split {
    input:
      inputs = attach_bcs_and_umis_main.output_bam,
      num_alignments = attach_bcs_and_umis_main.num_alignments
  }

  scatter(chunk_input in attach_bcs_and_umis_main.output_bam) {
    call bucket_by_bc_main {
      input:
        chunk_input = chunk_input,
        read_groups = bucket_by_bc_split.read_groups
    }
  }
  
  call sort_by_bc_main as sort_AA { input: bucketed_bams = bucket_by_bc_main.AA_bam }
  call sort_by_bc_main as sort_AT { input: bucketed_bams = bucket_by_bc_main.AT_bam }
  call sort_by_bc_main as sort_AG { input: bucketed_bams = bucket_by_bc_main.AG_bam }
  call sort_by_bc_main as sort_AC { input: bucketed_bams = bucket_by_bc_main.AC_bam }
  call sort_by_bc_main as sort_TA { input: bucketed_bams = bucket_by_bc_main.TA_bam }
  call sort_by_bc_main as sort_TT { input: bucketed_bams = bucket_by_bc_main.TT_bam }
  call sort_by_bc_main as sort_TG { input: bucketed_bams = bucket_by_bc_main.TG_bam }
  call sort_by_bc_main as sort_TC { input: bucketed_bams = bucket_by_bc_main.TC_bam }
  call sort_by_bc_main as sort_GA { input: bucketed_bams = bucket_by_bc_main.GA_bam }
  call sort_by_bc_main as sort_GT { input: bucketed_bams = bucket_by_bc_main.GT_bam }
  call sort_by_bc_main as sort_GG { input: bucketed_bams = bucket_by_bc_main.GG_bam }
  call sort_by_bc_main as sort_GC { input: bucketed_bams = bucket_by_bc_main.GC_bam }
  call sort_by_bc_main as sort_CA { input: bucketed_bams = bucket_by_bc_main.CA_bam }
  call sort_by_bc_main as sort_CT { input: bucketed_bams = bucket_by_bc_main.CT_bam }
  call sort_by_bc_main as sort_CG { input: bucketed_bams = bucket_by_bc_main.CG_bam }
  call sort_by_bc_main as sort_CC { input: bucketed_bams = bucket_by_bc_main.CC_bam }
  call sort_by_bc_main as sort_null { input: bucketed_bams = bucket_by_bc_main.null_bam }
  
  call sort_by_bc_join {
    input:
      sorted_AA = sort_AA.default, sorted_AT = sort_AT.default, sorted_AC = sort_AC.default, sorted_AG = sort_AG.default,
      sorted_TA = sort_TA.default, sorted_TT = sort_TT.default, sorted_TC = sort_TC.default, sorted_TG = sort_TG.default,
      sorted_CA = sort_CA.default, sorted_CT = sort_CT.default, sorted_CC = sort_CC.default, sorted_CG = sort_CG.default,
      sorted_GA = sort_GA.default, sorted_GT = sort_GT.default, sorted_GC = sort_GC.default, sorted_GG = sort_GG.default,
      sorted_null = sort_null.default,
      total_AA = sort_AA.total_reads, total_AT = sort_AT.total_reads, total_AC = sort_AC.total_reads, total_AG = sort_AG.total_reads,
      total_TA = sort_TA.total_reads, total_TT = sort_TT.total_reads, total_TC = sort_TC.total_reads, total_TG = sort_TG.total_reads,
      total_CA = sort_CA.total_reads, total_CT = sort_CT.total_reads, total_CC = sort_CC.total_reads, total_CG = sort_CG.total_reads,
      total_GA = sort_GA.total_reads, total_GT = sort_GT.total_reads, total_GC = sort_GC.total_reads, total_GG = sort_GG.total_reads,
      total_null = sort_null.total_reads
  }

  call mark_duplicates_split {
    input:
      input_bam = sort_by_bc_join.default
  }
  
  scatter(chunk in mark_duplicates_split.chunks) {
    call mark_duplicates_main {
      input:
        input_bam = sort_by_bc_join.default,
        chunk = chunk,
        align = extract_reads_join.out_align,
        reference_path = reference_path
    }
  }

  call mark_duplicates_join {
    input:
      output_bam_from_mains = mark_duplicates_main.output_bam,
      chunked_reporter_from_mains = mark_duplicates_main.chunked_reporter,
      outs = mark_duplicates_main.out
  }

  call count_genes_split {
    input:
      bam_inputs = mark_duplicates_main.output_bam,
      barcode_whitelist = setup_chunks.barcode_whitelist,
      gem_groups = extract_reads_join.gem_groups,
      reference_path = reference_path,
      barcode_summary = attach_bcs_and_umis_join.barcode_summary,
      chemistry_def = setup_chunks.chemistry_def
  }

  scatter(chunk in count_genes_split.chunks) {
    call count_genes_main {
    input:
      bam_inputs = mark_duplicates_main.output_bam,
      barcode_whitelist = setup_chunks.barcode_whitelist,
      gem_groups = extract_reads_join.gem_groups,
      reference_path = reference_path,
      barcode_summary = attach_bcs_and_umis_join.barcode_summary,
      chemistry_def = setup_chunks.chemistry_def,
      chunk = chunk,
      align = extract_reads_join.out_align
    }
  }

  call count_genes_join {
    input:
      matrices_h5_from_mains = count_genes_main.matrices_h5,
      chunked_reporter_from_mains = count_genes_main.chunked_reporter,
      chemistry_def = setup_chunks.chemistry_def,
      gem_groups = extract_reads_join.gem_groups,
      sample_id = sample_id,
      outs = count_genes_main.outs
  }
  
  call filter_barcodes {
    input:
      sample_id = sample_id,
      matrices_h5 = count_genes_join.matrices_h5,
      barcode_whitelist = setup_chunks.barcode_whitelist,
      gem_groups = extract_reads_join.gem_groups,
      chemistry_def = setup_chunks.chemistry_def,
      raw_fastq_summary = extract_reads_join.summary,
      attach_bcs_summary = attach_bcs_and_umis_join.summary,
      barcode_summary = attach_bcs_and_umis_join.barcode_summary
  }

  output {
    # summarize reports was not staged, do not have:
    # web_summary.html
    # metrics_summary.csv

    # some metrics we do have but aren't reported by cellranger:
    File attach_bcs_and_umis_summary=attach_bcs_and_umis_join.summary
    File filter_barcodes_summary=filter_barcodes.summary
    File count_genes_summary=count_genes_join.reporter_summary
    File extract_reads_summary=extract_reads_join.summary
    File mark_duplicates_summary=mark_duplicates_join.summary
    File raw_gene_bc_matrices_mex=count_genes_join.matrices_mex
    File raw_gene_bc_matrices_h5=count_genes_join.matrices_h5
    File filtered_gene_bc_matrices_mex=filter_barcodes.filtered_matrices_mex
    File filtered_gene_bc_matrices_h5=filter_barcodes.filtered_matrices_h5
    File bam_output=sort_by_bc_join.default
  }
}
