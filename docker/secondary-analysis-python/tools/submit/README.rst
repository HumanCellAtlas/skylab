Overview
========
This package provides utilities for submitting an analysis bundle to the Human Cell Atlas Data Coordination Platform.

The steps in the submission process are as follows:

* Create analysis.json
* Create submission envelope and upload metadata
* Get URN need to stage files
* Stage files
* Confirm submission

analysis_json.py
================
Creates analysis.json file.

Invoke it like this::

    analysis-json \
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

All arguments are required.

submit.py
=========
Creates submission envelope and uploads metadata.

Invoke it like this::

    submit \  
      -submit_url ${submit_url} \  
      -analysis_json_path analysis.json  

Both arguments are required.

submission_urn.py
=================
Obtains URN needed for staging files. Queries ingest API until URN is available.

Invoke it like this::

    submission-urn \
      -envelope_url ${submission_url} \
      -retry_seconds ${retry_seconds} \
      -timeout_seconds ${timeout_seconds} > submission_urn.txt

All arguments are required.

stage.py
========
Uploads files to staging area.

Invoke it like this::

    stage $local_file_path $submission_urn

Both arguments are required

confirm.py
==========
Confirms submission. This causes the ingest service to finalize the submission and create a bundle in the storage service.

Waits until submission status is "Valid", since submission cannot be confirmed until then.

Invoke it like this::

    confirm \
      -envelope_url ${submission_url} \
      -retry_seconds ${retry_seconds} \
      -timeout_seconds ${timeout_seconds}

All arguments are required
