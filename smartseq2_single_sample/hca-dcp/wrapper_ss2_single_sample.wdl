import "ss2_single_sample.wdl" as ss2
import "submit.wdl" as submit_wdl


task GetInputs {
  String bundle_uuid
  String bundle_version
  String dss_url
  Int retry_seconds
  Int timeout_seconds

  command <<<
    python <<CODE
    from secondary_analysis import utils

    # Get bundle manifest
    uuid = '${bundle_uuid}'
    version = '${bundle_version}'
    dss_url = '${dss_url}'
    retry_seconds = ${retry_seconds}
    timeout_seconds = ${timeout_seconds}
    print('Getting bundle manifest for id {0}, version {1}'.format(uuid, version))
    manifest_files = utils.get_manifest_files(uuid, version, dss_url, timeout_seconds, retry_seconds)

    print('Downloading assay.json')
    assay_json_uuid = manifest_files['name_to_meta']['assay.json']['uuid']
    assay_json = utils.get_file_by_uuid(assay_json_uuid, dss_url)

    sample_id = assay_json['sample_id']
    fastq_1_name = assay_json['seq']['lanes'][0]['r1']
    fastq_2_name = assay_json['seq']['lanes'][0]['r2']
    fastq_1_url = manifest_files['name_to_meta'][fastq_1_name]['url']
    fastq_2_url = manifest_files['name_to_meta'][fastq_2_name]['url']

    print('Creating input map')
    with open('inputs.tsv', 'w') as f:
        f.write('fastq_1\tfastq_2\tsample_id\n')
        f.write('{0}\t{1}\t{2}\n'.format(fastq_1_url, fastq_2_url, sample_id))
    print('Wrote input map')
    CODE
  >>>
  runtime {
    docker: "humancellatlas/secondary-analysis-python:0.1.1"
  }
  output {
    Object inputs = read_object("inputs.tsv")
  }
}

workflow WrapperSs2RsemSingleSample {
  String bundle_uuid
  String bundle_version

  File gtf
  File ref_fasta
  File rrna_interval
  File ref_flat
  String star_genome
  String rsem_genome
  String reference_bundle

  # Submission
  File format_map
  String dss_url
  String submit_url
  String method
  String schema_version
  String run_type
  Int retry_seconds
  Int timeout_seconds

  call GetInputs as prep {
    input:
      bundle_uuid = bundle_uuid,
      bundle_version = bundle_version,
      dss_url = dss_url,
      retry_seconds = retry_seconds,
      timeout_seconds = timeout_seconds
  }

  call ss2.Ss2RsemSingleSample as analysis {
    input:
      fastq_read1 = prep.inputs.fastq_1,
      fastq_read2 = prep.inputs.fastq_2,
      output_prefix = prep.inputs.sample_id,
      gtf = gtf,
      ref_fasta = ref_fasta,
      rrna_interval = rrna_interval,
      ref_flat = ref_flat,
      star_genome = star_genome,
      rsem_genome = rsem_genome
  }

  call submit_wdl.submit {
    input:
      inputs = [
        {
          'name': 'fastq_read1',
          'value': prep.inputs.fastq_1,
        },
        {
          'name': 'fastq_read2',
          'value': prep.inputs.fastq_2,
        },
        {
          'name': 'output_prefix',
          'value': prep.inputs.sample_id,
        },
        {
          'name': 'gtf',
          'value': gtf,
        },
        {
          'name': 'ref_fasta',
          'value': ref_fasta,
        },
        {
          'name': 'rrna_interval',
          'value': rrna_interval,
        },
        {
          'name': 'ref_flat',
          'value': ref_flat,
        },
        {
          'name': 'star_genome',
          'value': star_genome,
        },
        {
          'name': 'rsem_genome',
          'value': rsem_genome,
        }
      ],
      outputs = [
        analysis.bam_file,
        analysis.bam_trans,
        analysis.rna_metrics,
        analysis.aln_metrics,
        analysis.rsem_gene_results,
        analysis.rsem_isoform_results,
        analysis.rsem_gene_count,
        analysis.gene_unique_counts,
        analysis.exon_unique_counts,
        analysis.transcript_unique_counts
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
