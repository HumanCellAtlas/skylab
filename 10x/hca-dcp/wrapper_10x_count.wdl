import "count.wdl" as countwdl
import "submit_10x_count.wdl" as submit

task GetInputs {
  String bundle_uuid
  String bundle_version

  command <<<
    python <<CODE
    import json
    import requests
    import subprocess

    # Get bundle manifest
    uuid = '${bundle_uuid}'
    version = '${bundle_version}'
    print('Getting bundle manifest for id {0}, version {1}'.format(uuid, version))

    url = "https://dss.dev.data.humancellatlas.org/v1/bundles/" + uuid + "?version=" + version + "&replica=gcp&directurls=true"
    print('GET {0}'.format(url))
    response = requests.get(url)
    print('{0}'.format(response.status_code))
    print('{0}'.format(response.text))
    manifest = response.json()

    bundle = manifest['bundle']
    uuid_to_url = {}
    name_to_meta = {}
    for f in bundle['files']:
        uuid_to_url[f['uuid']] = f['url']
        name_to_meta[f['name']] = f

    print('Downloading assay.json')
    assay_json_uuid = name_to_meta['assay.json']['uuid']
    url = "https://dss.dev.data.humancellatlas.org/v1/files/" + assay_json_uuid + "?replica=gcp"
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

    print('Creating input map')
    with open('inputs.tsv', 'w') as f:
        f.write('r1\tr2\ti1\tsample_id\n')
        f.write('{0}\t{1}\t{2}\t{3}\n'.format(r1, r2, i1, sample_id))
    print('Wrote input map')
    CODE
  >>>
  runtime {
    docker: "humancellatlas/secondary-analysis-python"
  }
  output {
    Object inputs = read_object("inputs.tsv")
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

  call GetInputs as prep {
    input:
      bundle_uuid = bundle_uuid,
      bundle_version = bundle_version
  }

  call countwdl.count as analysis {
    input:
      sample_def = sample_def,
      r1 = prep.inputs.r1,
      r2 = prep.inputs.r2,
      i1 = prep.inputs.i1,
      sample_id = prep.inputs.sample_id,
      reads_per_file = reads_per_file,
      subsample_rate = subsample_rate,
      primers = primers,
      align = align,
      reference_path = reference_path,
      umi_min_qual_threshold = umi_min_qual_threshold
  }

  call submit.SubmitWrapper10xCount {
    input:
      attach_bcs_and_umis_summary = analysis.attach_bcs_and_umis_summary,
      filter_barcodes_summary = analysis.filter_barcodes.summary,
      count_genes_summary = analysis.count_genes_join.reporter_summary,
      extract_reads_summary = analysis.extract_reads_join.summary,
      mark_duplicates_summary = analysis.mark_duplicates_join.summary,
      raw_gene_bc_matrices_mex = analysis.count_genes_join.matrices_mex,
      raw_gene_bc_matrices_h5 = analysis.count_genes_join.matrices_h5,
      filtered_gene_bc_matrices_mex = analysis.filter_barcodes.filtered_matrices_mex,
      filtered_gene_bc_matrices_h5 = analysis.filter_barcodes.filtered_matrices_h5,
      bam_output = analysis.sort_by_bc_join.default

  }
}
