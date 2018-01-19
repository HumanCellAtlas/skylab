# The two preflight stages mostly seem to do some basic sanity checks: verification of
# parameters, things like that.
task make_fastqs_preflight_local {
  File run_path
  Boolean check_executables
  Array[Int]? lanes
  Array[Map[String, File]] specs
  String barcode_whitelist
  String bc_read_type
  Int bc_start_index
  Int bc_length
  String si_read_type
  String umi_read_type
  Int umi_start_index
  Int umi_length
  String bcl2fastq2_args
  String? bases_mask
  Boolean ignore_dual_index

  command {
    /opt/tenx/run_in_10x_env.bash martian stage run MAKE_FASTQS_PREFLIGHT_LOCAL main \
      --run_path '${run_path}' \
      --check_executables ${true='true' false='false' check_executables} \
      ${'--lanes ' + lanes} \
      --specs '${sep="\' \'" specs}' \
      --barcode_whitelist '${barcode_whitelist}' \
      --bc_read_type '${bc_read_type}' \
      --bc_start_index '${bc_start_index}' \
      --bc_length '${bc_length}' \
      --si_read_type '${si_read_type}' \
      --umi_read_type '${umi_read_type}' \
      --umi_length '${umi_length}' \
      --bcl2fastq2_args '${bcl2fastq2_args}' \
      ${'--bases_mask ' + bases_mask} \
      --ignore_dual_index ${true='true' false='false' check_executables}
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File make_fastqs_preflight_local_stdout = stdout()
  }
}

task make_fastqs_preflight {
  File run_path
  File? output_path
  File? interop_output_path
  String barcode_whitelist
  Boolean check_executables
  Int max_bcl2fastq_threads

  command {
    /opt/tenx/run_in_10x_env.bash martian stage run MAKE_FASTQS_PREFLIGHT main \
      --run_path '${run_path}' \
      --check_executables ${true='true' false='false' check_executables} \
      ${'--output_path ' +  output_path} \
      ${'--interop_output_path ' +  interop_output_path} \
      --barcode_whitelist '${barcode_whitelist}' \
      --max_bcl2fastq_threads '${max_bcl2fastq_threads}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File make_fastqs_preflight_local_stdout = stdout()
  }

}

# Prepare samplesheet creates the special 10x samplesheet for bcl2fastq
task prepare_samplesheet {
  File run_path
  Array[Map[String, File]] specs
  String project
  String bc_read_type
  Int bc_length
  String si_read_type

  command {
    /opt/tenx/run_in_10x_env.bash martian stage run PREPARE_SAMPLESHEET main \
      --run_path '${run_path}' \
      --bc_read_type '${bc_read_type}' \
      --project '${project}' \
      --bc_length '${bc_length}' \
      --si_read_type '${si_read_type}' \
      --specs '${sep="\' \' " specs}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }
  
  output {
    File samplesheet = "samplesheet.csv"
    File input_samplesheet = "input_samplesheet.csv"
    Boolean dual_indexed_samplesheet = read_boolean("dual_indexed_samplesheet")
  }

}

# This runs bcl2fastq with the samplesheet created above
task bcl2fastq_with_samplesheet_main {
  File run_path
  File? output_path
  File? interop_output_path
  File samplesheet_path
  String? bases_mask
  String si_read_type
  String? bcl2fastq1_args
  String bcl2fastq2_args
  Int max_bcl2fastq_threads
  Boolean dual_indexed_samplesheet
  Boolean ignore_dual_index
  
  command {
    /opt/tenx/run_in_10x_env.bash martian stage run BCL2FASTQ_WITH_SAMPLESHEET main \
      --run_path '${run_path}' \
      ${'--output_path ' +  output_path} \
      ${'--interop_output_path ' +  interop_output_path} \
      --samplesheet_path '${samplesheet_path}' \
      ${'--bases_mask ' + bases_mask} \
      --si_read_type '${si_read_type}' \
      ${'--bcl2fastq1_args ' + bcl2fastq1_args} \
      --bcl2fastq2_args '${bcl2fastq2_args}' \
      --max_bcl2fastq_threads '${max_bcl2fastq_threads}' \
      --dual_indexed_samplesheet ${true='true' false='false' dual_indexed_samplesheet} \
      --ignore_dual_index ${true='true' false='false' ignore_dual_index}
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Boolean rc_i2_read = read_boolean("rc_i2_read")
    Map[String, String] file_read_types_map = read_map("file_read_types_map")
    String bcl2fastq_version = read_string("bcl2fastq_version")
    String bcl2fastq_args = read_string("bcl2fastq_args")
    File fastq_path = "fastq_path"
    File interop_path = "interop_path"
  }
}

# This is a little weird since it's a join without a split. But the split in the 10x code is
# basically a no-op, so I've dropped it from the WDL. Moreover, this join isn't really
# joining anything, it's more of a continuation of main.
task bcl2fastq_with_samplesheet_join {
  File run_path
  File? output_path
  File? interop_output_path
  File samplesheet_path
  String? bases_mask
  String si_read_type
  String? bcl2fastq1_args
  String bcl2fastq2_args
  Int max_bcl2fastq_threads
  Boolean dual_indexed_samplesheet
  Boolean ignore_dual_index

  File fastq_path
  File interop_path
  Boolean rc_i2_read_chunk
  Map[String, String] file_read_types_map_chunk
  String bcl2fastq_version_chunk
  String bcl2fastq_args_chunk
  
  command {
    /opt/tenx/run_in_10x_env.bash martian stage run BCL2FASTQ_WITH_SAMPLESHEET join \
      --run_path '${run_path}' \
      ${'--output_path ' +  output_path} \
      ${'--interop_output_path ' +  interop_output_path} \
      --samplesheet_path '${samplesheet_path}' \
      ${'--bases_mask ' + bases_mask} \
      --si_read_type '${si_read_type}' \
      --bcl2fastq1_args '${bcl2fastq1_args}' \
      --bcl2fastq2_args '${bcl2fastq2_args}' \
      --max_bcl2fastq_threads '${max_bcl2fastq_threads}' \
      --dual_indexed_samplesheet ${true='true' false='false' dual_indexed_samplesheet} \
      --ignore_dual_index ${true='true' false='false' ignore_dual_index} \
      --fastq_path_output '${fastq_path}' \
      --interop_path_output '${interop_path}' \
      --rc_i2_read_output  ${true='true' false='false' rc_i2_read_chunk} \
      --file_read_types_map_output ${write_map(file_read_types_map_chunk)} \
      --bcl2fastq_version_output '${bcl2fastq_version_chunk}' \
      --bcl2fastq_args_output '${bcl2fastq_args_chunk}'

  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Boolean rc_i2_read = read_boolean("rc_i2_read")
    Map[String, String] file_read_types_map = read_map("file_read_types_map")
    String bcl2fastq_version = read_string("bcl2fastq_version")
    String bcl2fastq_args = read_string("bcl2fastq_args")
  }
}

# This is the first real split/main/join step. It creates some JSON files with
# quality metrics.
task make_qc_summary_split {
  File fastq_path
  
  # Cell Ranger wants to work with paths that persists from call to calla lot, but the
  # /cromwell-execution/etc/etc path we get from File inputs don't do that. So hard
  # link the path to $CWD/fastq_path before we start working with it.
  command {
    cp -al '${fastq_path}' fastq_path

    /opt/tenx/run_in_10x_env.bash martian stage run MAKE_QC_SUMMARY split \
      --fastq_path fastq_path

    for split in martian_split_*; do
      jq -r '.input_files | join("\t")' "$split" >> split_input_files
    done
  }
  
  runtime {
    docker: "10x_hca_wdl:latest"
  }

  # The output of the split step is a list of dicts. Each dict has a keys that are
  # strings, and values of various types. In this case, the values are strings and
  # ints. There's also an "input_files" field that martian says is of type string,
  # but that's a trick, a ruse. It's actually a list of files, and we have to
  # break it out separately or else Cromwell won't know to output the files and
  # stage them for future steps.
  output {
    Array[File] split_files = glob("martian_split_*")
    Array[Array[File]] input_files = read_tsv("split_input_files")
  }
}

task make_qc_summary_main {
  File run_path
  File fastq_path
  File interop_path
  String barcode_whitelist
  String bc_read_type
  Int bc_start_index
  Int bc_length
  String si_read_type
  String umi_read_type
  Int umi_start_index
  Int umi_length
  Boolean rc_i2_read
  Map[String, String] file_read_types_map
  String software_version
  String bcl2fastq_version
  String bcl2fastq_args
  File split_file
  Array[File] input_files
  
  # If you want to use bash tricks in a WDL command, this is apparently what you have
  # to do. Otherwise WDL will try to interpret it and complain.
  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"
  
  # Note again that the paths we get from input_files are going to confuse Cell Ranger,
  # so link them to $CWD/fastq_path
  command {
    mkdir fastq_path
    for fastq in '${sep="\' \'" input_files}'; do
      mkdir -p $(dirname fastq_path/${dollar}${lbrace}fastq##*/fastq_path${rbrace})
      ln $fastq fastq_path/${dollar}${lbrace}fastq##*/fastq_path${rbrace}
    done

    /opt/tenx/run_in_10x_env.bash martian stage run MAKE_QC_SUMMARY main \
      --run_path '${run_path}' \
      --fastq_path fastq_path \
      --interop_path '${interop_path}' \
      --barcode_whitelist '${barcode_whitelist}' \
      --bc_read_type '${bc_read_type}' \
      --bc_start_index '${bc_start_index}' \
      --bc_length '${bc_length}' \
      --si_read_type '${si_read_type}' \
      --umi_read_type '${umi_read_type}' \
      --umi_length '${umi_length}' \
      --rc_i2_read  ${true='true' false='false' rc_i2_read} \
      --file_read_types_map ${write_map(file_read_types_map)} \
      --software_version '${software_version}' \
      --bcl2fastq_version '${bcl2fastq_version}' \
      --bcl2fastq_args '${bcl2fastq_args}' \
      --split_file '${split_file}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }
  
  # The MRO for this stage says it should output a json file called qc_summary. But that's just
  # another sneaky trick! It actually outputs a json string that points to a bunch of local json
  # files. So just output those files so Cromwell knows we need them for later.
  output {
    Array[File] qc_jsons = glob("output_*.json")
  }
}

task make_qc_summary_join {
  
  File run_path
  String barcode_whitelist
  String bc_read_type
  Int bc_start_index
  Int bc_length
  String si_read_type
  String umi_read_type
  Int umi_start_index
  Int umi_length
  String software_version
  String bcl2fastq_version
  String bcl2fastq_args

  Array[Array[File]] qc_jsons
  Array[File] split_files
  
  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"

  command {

    samples="${dollar}(jq -s '.[].sample' '${sep="\' \'" split_files}' | jq -rs 'join(" ")')"
    lanes="${dollar}(jq -s '.[].lane | tostring' '${sep="\' \'" split_files} | jq -rs 'join(" ")')"

    # The mro for this stage says that the output of the main steps should be a json file called
    # qc summary, so that's what the martian_cli is looking for for the join step. But the main
    # steps don't actually output that, so we have to make our own here. This code iterates
    # through the collections of files that were produced by the main steps and creates json
    # files with the expected structure. The json should have three keys: barcode, read1, and
    # read2. And the values for each is a list of paths. You can tell which list each path should
    # go it by the presence of a substring: BC --> barcode, R1 --> read1, R2 --> read2
    echo '[${sep=',' qc_jsons}]' | \
      jq -c '.[] | ${lbrace}"barcode": map(select( . |contains("BC"))), "read1": map(select( . |contains("R1"))), "read2": map(select( . |contains("R2")))${rbrace}' | \
      split -l 1 -a 2 - qc_json_split

    echo "${dollar}${lbrace}lanes${rbrace}"

    /opt/tenx/run_in_10x_env.bash martian stage run MAKE_QC_SUMMARY join \
      --run_path '${run_path}' \
      --sample_split ${dollar}${lbrace}samples${rbrace} \
      --lane_split ${dollar}${lbrace}lanes${rbrace} \
      --qc_summary_output qc_json_split* \
      --barcode_whitelist '${barcode_whitelist}' \
      --bc_read_type '${bc_read_type}' \
      --bc_start_index '${bc_start_index}' \
      --bc_length '${bc_length}' \
      --si_read_type '${si_read_type}' \
      --umi_read_type '${umi_read_type}' \
      --umi_start_index '${umi_start_index}' \
      --umi_length '${umi_length}' \
      --software_version '${software_version}' \
      --bcl2fastq_version '${bcl2fastq_version}' \
      --bcl2fastq_args '${bcl2fastq_args}'
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    File qc_summary = "qc_summary.json"
  }
}

task merge_fastqs_by_lane_sample_split {
  
  File fastq_path
  File samplesheet_path
  String bcl2fastq_version

  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"

  command {
    cp -al '${fastq_path}' fastq_path
    /opt/tenx/run_in_10x_env.bash martian stage run MERGE_FASTQS_BY_LANE_SAMPLE split \
      --fastq_path fastq_path \
      --samplesheet_path '${samplesheet_path}' \
      --bcl2fastq_version '${bcl2fastq_version}'

    for split in martian_split_*; do
      jq -r '.input_files | join("\t")' "$split" >> split_input_files
    done
  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] split_files = glob("martian_split_*")
    Array[Array[File]] input_files = read_tsv("split_input_files")
  }
}

task merge_fastqs_by_lane_sample_main {

  File samplesheet_path
  String bcl2fastq_version
  File split_file
  Array[File] input_files

  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"

  command {
    # Keep the fastq paths as $PWD/fastq_path so when the paths are passed to the join
    # step, they'll still be valid.
    mkdir fastq_path
    for fastq in '${sep="\' \'" input_files}'; do
      mkdir -p $(dirname fastq_path/${dollar}${lbrace}fastq##*/fastq_path${rbrace})
      ln $fastq fastq_path/${dollar}${lbrace}fastq##*/fastq_path${rbrace}
    done

    /opt/tenx/run_in_10x_env.bash martian stage run MERGE_FASTQS_BY_LANE_SAMPLE main \
      --fastq_path fastq_path \
      --samplesheet_path '${samplesheet_path}' \
      --bcl2fastq_version '${bcl2fastq_version}' \
      --split_file '${split_file}'

    # The output of this step is "merged_file_paths", which in martian is a string[],
    # but in WDL, we really want an Array[File]. The martian_cli writes a json file
    # called merged_file_paths that's a list of paths. Cromwell can't read that, since
    # read_json is not yet implemented. So this line converts the json to a
    # newline-separated list of paths that works with read_lines.
    jq -r '.[]' merged_file_paths > merged_file_paths.lines
  }
  
  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] merged_file_paths = read_lines("merged_file_paths.lines")
    Boolean files_merged = read_boolean("files_merged")
  }
}

task merge_fastqs_by_lane_sample_join {

  Array[Array[File]] merged_file_paths_output

  String dollar = "$"
  String lbrace = "{"
  String rbrace = "}"

  command {
    
    # merged_file_paths_output is an Array[Array[File]], so we have to refer
    # to it with a "sep" or cromwell will complain. When we sep with "," and
    # stick a couple brackets around it, we get a string that is a json list
    # of lists of paths. We then iterate over each top level list and for each
    # path in each sublist, we remove the cromwell-specific part of the path by
    # deleting everything up to "fastq_path". Later, we will create links in
    # $PWD/fastq_path so these truncated paths point to something valid.
    #
    # After truncating the paths, we write each list of paths to its own file.
    # The `jq -c` will print each list as a single line, and the `split -l 1`
    # will print each line to its own file, prefixed by "merged_file_paths_output"
    echo '[${sep=',' merged_file_paths_output}]' | \
      jq -c '.[] | map(sub(".*(?=fastq_path\/)"; "")) ' | \
      split -l 1 -a 2 - merged_file_paths_output
    
    # Now we are once again iterating over the merged_file_paths_output variable.
    # The `echo | jq .[] | .[]` will pass each file path into this loop as the
    # "fastq" variable. Then, for each path, strip off the cromwell-specific
    # path of the path, just like above, and link to the $PWD/fastq_path
    # directory
    while read fastq; do
      echo ${dollar}${lbrace}fastq${rbrace}
      mkdir -p $(dirname fastq_path/${dollar}${lbrace}fastq##*/fastq_path${rbrace})
      ln $fastq fastq_path/${dollar}${lbrace}fastq##*/fastq_path${rbrace}
    done < <( echo '[${sep=',' merged_file_paths_output}]' | jq -r '.[] | .[]')


    /opt/tenx/run_in_10x_env.bash martian stage run MERGE_FASTQS_BY_LANE_SAMPLE join \
      --merged_file_paths_output merged_file_paths_output*

    # The join step just renames files in fastq_path, so we have to find the outputs
    # ourselves. This can be done by just finding all the files in fastq_path and
    # writing them to a file that read_lines and read.
    find fastq_path -type f -exec echo ${lbrace}${rbrace} >> merged_file_paths.lines \;

  }

  runtime {
    docker: "10x_hca_wdl:latest"
  }

  output {
    Array[File] merged_file_paths = read_lines("merged_file_paths.lines")
  }
}

workflow make_fastqs {
  File run_path
  Boolean check_executables
  Array[Int]? lanes
  Array[Map[String, File]] specs
  String barcode_whitelist
  String bc_read_type
  Int bc_start_index
  Int bc_length
  String si_read_type
  String umi_read_type
  Int umi_start_index
  Int umi_length
  String? bcl2fastq1_args
  String bcl2fastq2_args
  String? bases_mask
  Boolean ignore_dual_index
  File? output_path
  File? interop_output_path
  Int max_bcl2fastq_threads
  String project
  String software_version

  call make_fastqs_preflight_local {
    input: run_path=run_path, check_executables=false, lanes=lanes, specs=specs, barcode_whitelist=barcode_whitelist,
      bc_read_type=bc_read_type, bc_start_index=bc_start_index, bc_length=bc_length, si_read_type=si_read_type,
      umi_read_type=umi_read_type, umi_start_index=umi_start_index, umi_length=umi_length, bcl2fastq2_args=bcl2fastq2_args,
      bases_mask=bases_mask, ignore_dual_index=ignore_dual_index
  }
  call make_fastqs_preflight {
    input: run_path=run_path, output_path=output_path, interop_output_path=interop_output_path, barcode_whitelist=barcode_whitelist,
      check_executables=true, max_bcl2fastq_threads=max_bcl2fastq_threads
  }

  call prepare_samplesheet {
    input: run_path=run_path, specs=specs, project=project, bc_read_type=bc_read_type, bc_length=bc_length, si_read_type=si_read_type
  }

  call bcl2fastq_with_samplesheet_main {
    input: run_path=run_path, output_path=output_path, interop_output_path=interop_output_path, samplesheet_path=prepare_samplesheet.samplesheet,
      bases_mask=bases_mask, si_read_type=si_read_type, bcl2fastq1_args=bcl2fastq1_args, bcl2fastq2_args=bcl2fastq2_args,
      max_bcl2fastq_threads=max_bcl2fastq_threads, dual_indexed_samplesheet=prepare_samplesheet.dual_indexed_samplesheet,
      ignore_dual_index=ignore_dual_index
  }

  call bcl2fastq_with_samplesheet_join {
    input: run_path=run_path, output_path=output_path, interop_output_path=interop_output_path, samplesheet_path=prepare_samplesheet.samplesheet,
      bases_mask=bases_mask, si_read_type=si_read_type, bcl2fastq1_args=bcl2fastq1_args, bcl2fastq2_args=bcl2fastq2_args,
      max_bcl2fastq_threads=max_bcl2fastq_threads, dual_indexed_samplesheet=prepare_samplesheet.dual_indexed_samplesheet,
      ignore_dual_index=ignore_dual_index, fastq_path=bcl2fastq_with_samplesheet_main.fastq_path,
      interop_path=bcl2fastq_with_samplesheet_main.interop_path, rc_i2_read_chunk=bcl2fastq_with_samplesheet_main.rc_i2_read,
      file_read_types_map_chunk=bcl2fastq_with_samplesheet_main.file_read_types_map,
      bcl2fastq_version_chunk=bcl2fastq_with_samplesheet_main.bcl2fastq_version,
      bcl2fastq_args_chunk=bcl2fastq_with_samplesheet_main.bcl2fastq_args
  }

  call make_qc_summary_split {
    input: fastq_path=bcl2fastq_with_samplesheet_main.fastq_path
  }

  scatter(split_pair in zip(make_qc_summary_split.split_files, make_qc_summary_split.input_files)) {
    call make_qc_summary_main {
      input: run_path=run_path, fastq_path=bcl2fastq_with_samplesheet_main.fastq_path, interop_path=bcl2fastq_with_samplesheet_main.interop_path,
        barcode_whitelist=barcode_whitelist, bc_read_type=bc_read_type, bc_length=bc_length, si_read_type=si_read_type, umi_read_type=umi_read_type,
        umi_length=umi_length, rc_i2_read=bcl2fastq_with_samplesheet_join.rc_i2_read, file_read_types_map=bcl2fastq_with_samplesheet_join.file_read_types_map,
        software_version=software_version, bcl2fastq_version=bcl2fastq_with_samplesheet_join.bcl2fastq_version,
        bcl2fastq_args=bcl2fastq_with_samplesheet_join.bcl2fastq_args, umi_start_index=umi_start_index, bc_start_index=bc_start_index,
        split_file=split_pair.left, input_files=split_pair.right
    }
  }

  call make_qc_summary_join {
    input: run_path=run_path, barcode_whitelist=barcode_whitelist, bc_read_type=bc_read_type, bc_start_index=bc_start_index, bc_length=bc_length,
      si_read_type=si_read_type, umi_read_type=umi_read_type, umi_start_index=umi_start_index, umi_length=umi_length,
      software_version=software_version, bcl2fastq_version=bcl2fastq_with_samplesheet_main.bcl2fastq_version,
      bcl2fastq_args=bcl2fastq_with_samplesheet_main.bcl2fastq_args,
      qc_jsons=make_qc_summary_main.qc_jsons, split_files=make_qc_summary_split.split_files
  }

  call merge_fastqs_by_lane_sample_split {
    input: fastq_path=bcl2fastq_with_samplesheet_main.fastq_path, samplesheet_path=prepare_samplesheet.samplesheet,
    bcl2fastq_version=bcl2fastq_with_samplesheet_main.bcl2fastq_version
  }

  scatter(split_pair in zip(merge_fastqs_by_lane_sample_split.split_files, merge_fastqs_by_lane_sample_split.input_files)) {
    call merge_fastqs_by_lane_sample_main {
      input: samplesheet_path=prepare_samplesheet.samplesheet, bcl2fastq_version=bcl2fastq_with_samplesheet_main.bcl2fastq_version,
        split_file=split_pair.left, input_files=split_pair.right
    }
  }

  call merge_fastqs_by_lane_sample_join {
    input: merged_file_paths_output=merge_fastqs_by_lane_sample_main.merged_file_paths
  }
}
