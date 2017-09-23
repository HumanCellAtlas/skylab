task GetMetadata {
  String analysis_output_path

  command <<<

    # Get workflow_id
    python <<CODE > workflow_id.txt

    url = '${analysis_output_path}'
    hash_end = url.rfind("/call-")
    hash_start = url.rfind('/', 0, hash_end) + 1
    hash = url[hash_start:hash_end]
    print(hash)

    CODE

    # Get metadata
    creds=/analysis-json/cromwell_credentials.txt
    curl -u $(cut -f1 $creds):$(cut -f2 $creds) \
      --compressed \
      "$(cut -f3 $creds)/api/workflows/v1/$(cat workflow_id.txt)/metadata" > metadata.json
  >>>
  runtime {
    docker: "gcr.io/broad-dsde-mint-dev/analysis-json:latest"
  }
  output {
    File metadata = "metadata.json"
    String workflow_id = read_string("workflow_id.txt")
  }
}

task create_submission {
  String workflow_id
  File metadata_json
  String bundle_uuid
  String reference_bundle
  String run_type
  String method
  String schema_version
  Array[Object] inputs
  Array[String] outputs
  File format_map
  String submit_url

  command <<<
    create-analysis-json \
      -analysis_id ${workflow_id} \
      -metadata_json ${metadata_json} \
      -input_bundles ${bundle_uuid} \
      -reference_bundle ${reference_bundle} \
      -run_type ${run_type} \
      -method ${method} \
      -schema_version ${schema_version} \
      -inputs_file ${write_objects(inputs)} \
      -outputs_file ${write_lines(outputs)} \
      -format_map ${format_map}

    create-envelope \
      -submit_url ${submit_url} \
      -analysis_json_path analysis.json
  >>>

  runtime {
    docker: "humancellatlas/secondary-analysis-python:test7"
  }
  output {
    File analysis_json = "analysis.json"
    String submission_url = read_string("submission_url.txt")
  }
}

task stage_and_confirm {
  String submission_url
  Array[File] files
  Int retry_seconds
  Int timeout_seconds
  String lb = "{"
  String rb = "}"

  command <<<
    staging_urn=$(get-staging-urn \
        -envelope_url ${submission_url} \
        -retry_seconds ${retry_seconds} \
        -timeout_seconds ${timeout_seconds})

    files=( ${sep=' ' files} )
    for f in "$${lb}files[@]${rb}"
    do
      echo "stage -d staging $f $staging_urn"
      stage -d staging $f $staging_urn
    done

    confirm-submission \
      -envelope_url ${submission_url} \
      -retry_seconds ${retry_seconds} \
      -timeout_seconds ${timeout_seconds}
  >>>

  runtime {
    docker: "humancellatlas/secondary-analysis-python:test7"
  }
}

workflow Submit {
  Array[Object] inputs
  Array[File] outputs
  File format_map
  String submit_url
  String bundle_uuid
  String reference_bundle
  String run_type
  String schema_version
  String method
  Int retry_seconds
  Int timeout_seconds

  call GetMetadata {
    input:
      analysis_output_path = outputs[0]
  }

  call create_submission {
    input:
      reference_bundle = reference_bundle,
      run_type = run_type,
      schema_version = schema_version,
      method = method,
      submit_url = submit_url,
      inputs = inputs,
      outputs = outputs,
      format_map = format_map,
      metadata_json = GetMetadata.metadata,
      bundle_uuid = bundle_uuid,
      workflow_id = GetMetadata.workflow_id,
  }

  call stage_and_confirm {
    input:
      submission_url = create_submission.submission_url,
      files = outputs,
      retry_seconds = retry_seconds,
      timeout_seconds = timeout_seconds
  }

  output {
    File analysis_json = create_submission.analysis_json
  }
}
