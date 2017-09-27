import "count.wdl" as countwdl
import "submit.wdl" as submit_wdl

task GetInputs {
  String bundle_uuid
  String bundle_version
  String dss_url
  Int retry_seconds
  Int timeout_seconds

  command <<<
    python <<CODE
    import json
    import requests
    import subprocess
    import time

    # Get bundle manifest
    uuid = '${bundle_uuid}'
    version = '${bundle_version}'
    print('Getting bundle manifest for id {0}, version {1}'.format(uuid, version))

    url = "${dss_url}/bundles/" + uuid + "?version=" + version + "&replica=gcp&directurls=true"
    start = time.time()
    current = start
    while current - start < ${timeout_seconds}:
        print('GET {0}'.format(url))
        response = requests.get(url)
        print('{0}'.format(response.status_code))
        print('{0}'.format(response.text))
        if 200 <= response.status_code <= 299:
            break
        time.sleep(${retry_seconds})
        current = time.time()
    manifest = response.json()

    bundle = manifest['bundle']
    uuid_to_url = {}
    name_to_meta = {}
    url_to_name_all = {}
    for f in bundle['files']:
        uuid_to_url[f['uuid']] = f['url']
        name_to_meta[f['name']] = f
        url_to_name_all[f['url']] = f['name']

    print('Downloading assay.json')
    assay_json_uuid = name_to_meta['assay.json']['uuid']
    url = "${dss_url}/files/" + assay_json_uuid + "?replica=gcp"
    print('GET {0}'.format(url))
    response = requests.get(url)
    print('{0}'.format(response.status_code))
    print('{0}'.format(response.text))
    assay_json = response.json()

    # Parse inputs from assay_json and write to inputs.tsv file
    sample_id = assay_json['sample_id']
    lanes = assay_json['seq']['lanes']
    r1 = [name_to_meta[lane['r1']]['url'] for lane in lanes]
    r2 = [name_to_meta[lane['r2']]['url'] for lane in lanes]
    i1 = [name_to_meta[lane['i1']]['url'] for lane in lanes]
    url_to_name = {}
    for i in r1:
        url_to_name[i] = url_to_name_all[i]
    for i in r2:
        url_to_name[i] = url_to_name_all[i]
    for i in i1:
        url_to_name[i] = url_to_name_all[i]

    with open('lanes.txt', 'w') as f:
        for i in range(len(lanes)):
            f.write('{}\n'.format(i))

    with open('r1.tsv', 'w') as f:
        for r in r1:
            f.write('{0}\n'.format(r))
    with open('r2.tsv', 'w') as f:
        for r in r2:
            f.write('{0}\n'.format(r))
    with open('i1.tsv', 'w') as f:
        for i in i1:
            f.write('{0}\n'.format(i))

    with open('r1_names.tsv', 'w') as f:
        for r in r1:
            f.write('{0}\n'.format(url_to_name[r]))
    with open('r2_names.tsv', 'w') as f:
        for r in r2:
            f.write('{0}\n'.format(url_to_name[r]))
    with open('i1_names.tsv', 'w') as f:
        for i in i1:
            f.write('{0}\n'.format(url_to_name[i]))
    print('Creating input map')
    with open('inputs.tsv', 'w') as f:
        f.write('sample_id\n')
        f.write('{0}\n'.format(sample_id))
    print('Wrote input map')
    CODE
  >>>
  runtime {
    docker: "humancellatlas/secondary-analysis-python"
  }
  output {
    Array[File] r1 = read_lines("r1.tsv")
    Array[File] r2 = read_lines("r2.tsv")
    Array[File] i1 = read_lines("i1.tsv")
    Array[String] r1_names = read_lines("r1_names.tsv")
    Array[String] r2_names = read_lines("r2_names.tsv")
    Array[String] i1_names = read_lines("i1_names.tsv")
    Array[Int] lanes = read_lines("lanes.txt")
    Object inputs = read_object("inputs.tsv")
  }
}

task rename_files {
  File r1
  File r2
  File i1
  String r1_name
  String r2_name
  String i1_name

  command <<<
    python <<CODE

    import subprocess

    subprocess.check_output(['mv', '${r1}', '${r1_name}'])
    subprocess.check_output(['mv', '${r2}', '${r2_name}'])
    subprocess.check_output(['mv', '${i1}', '${i1_name}'])

    CODE
  >>>
  runtime {
    docker: "humancellatlas/secondary-analysis-python"
  }
  output {
    File r1_new = "${r1_name}"
    File r2_new = "${r2_name}"
    File i1_new = "${i1_name}"
  }
}

task inputs_for_submit {
  Array[String] r1
  Array[String] r2
  Array[String] i1
  Array[Object] other_inputs
  Array[Object] primers

  command <<<
    primers_file="${write_objects(primers)}"
    mv $primers_file primers.tsv

    python <<CODE
    inputs = []
    print('other inputs')
    with open('${write_objects(other_inputs)}') as f:
        keys = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            input = {}
            for i, key in enumerate(keys):
                input[key] = values[i]
            print(input)
            inputs.append(input)

    print('r1')
    r1 = ['${sep="', '" r1}']
    for r in r1:
        input = {
            'name': r.split('\t')[-1],
            'value': r
        }
        inputs.append(input)

    print('r2')
    r2 = ['${sep="', '" r2}']
    for r in r2:
        input = {
            'name': r.split('\t')[-1],
            'value': r
        }
        inputs.append(input)

    print('i1')
    i1 = ['${sep="', '" i1}']
    for i in i1:
        input = {
            'name': i.split('\t')[-1],
            'value': i
        }
        inputs.append(input)

    print('write inputs.tsv')
    with open('inputs.tsv', 'w') as f:
        f.write('name\tvalue\n')
        for input in inputs:
            print(input)
            f.write('{0}\t{1}\n'.format(input['name'], input['value']))
    print('finished')
    CODE
  >>>
  runtime {
    docker: "humancellatlas/secondary-analysis-python"
  }
  output {
    File inputs = "inputs.tsv"
    File primers_tsv = "primers.tsv"
  }
}


workflow Wrapper10xCount {
  String bundle_uuid
  String bundle_version

  File sample_def
  Int reads_per_file
  Float subsample_rate
  Array[Map[String, String]] primers
  String align
  File reference_path
  Int umi_min_qual_threshold

  # Submission
  File format_map
  String dss_url
  String submit_url
  String reference_bundle
  String run_type
  String schema_version
  String method
  Int retry_seconds
  Int timeout_seconds

  call GetInputs {
    input:
      bundle_uuid = bundle_uuid,
      bundle_version = bundle_version,
      dss_url = dss_url,
      retry_seconds = retry_seconds,
      timeout_seconds = timeout_seconds
  }

  scatter(i in GetInputs.lanes) {
    call rename_files as prep {
      input:
        r1 = GetInputs.r1[i],
        r2 = GetInputs.r2[i],
        i1 = GetInputs.i1[i],
        r1_name = GetInputs.r1_names[i],
        r2_name = GetInputs.r2_names[i],
        i1_name = GetInputs.i1_names[i]
    }
  }

  call countwdl.count as analysis {
    input:
      sample_def = sample_def,
      r1 = prep.r1_new,
      r2 = prep.r2_new,
      i1 = prep.i1_new,
      sample_id = GetInputs.inputs.sample_id,
      reads_per_file = reads_per_file,
      subsample_rate = subsample_rate,
      primers = primers,
      align = align,
      reference_path = reference_path,
      umi_min_qual_threshold = umi_min_qual_threshold
  }

  String sample_id = GetInputs.inputs.sample_id
  File primers_tsv = inputs_for_submit.primers_tsv

  call inputs_for_submit {
    input:
      r1 = GetInputs.r1,
      r2 = GetInputs.r2,
      i1 = GetInputs.i1,
      primers = primers,
      other_inputs = [
        {
          'name': 'sample_def',
          'value': sample_def
        },
        {
          'name': 'sample_id',
          'value': sample_id
        },
        {
          'name': 'reads_per_file',
          'value': reads_per_file
        },
        {
          'name': 'subsample_rate',
          'value': subsample_rate
        },
        {
          'name': 'align',
          'value': align
        },
        {
          'name': 'reference_path',
          'value': reference_path
        },
        {
          'name': 'umi_min_qual_threshold',
          'value': umi_min_qual_threshold
        }
      ]
  }

  Array[Object] inputs = read_objects(inputs_for_submit.inputs)

  call submit_wdl.submit {
    input:
      inputs = inputs,
      outputs = [
#        analysis.attach_bcs_and_umis_summary,
#        analysis.filter_barcodes_summary,
        analysis.count_genes_summary,
#        analysis.extract_reads_summary,
#        analysis.mark_duplicates_summary,
        analysis.raw_gene_bc_matrices_mex,
        analysis.raw_gene_bc_matrices_h5,
        analysis.filtered_gene_bc_matrices_mex,
        analysis.filtered_gene_bc_matrices_h5,
        analysis.bam_output
      ],
      format_map = format_map,
      submit_url = submit_url,
      input_bundle_uuid = bundle_uuid,
      reference_bundle = reference_bundle,
      run_type = run_type,
      schema_version = schema_version,
      method = method,
      retry_seconds = retry_seconds,
      timeout_seconds = timeout_seconds
  }
}
