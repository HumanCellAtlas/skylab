import "ss2_single_sample.wdl" as ss2
import "submit_ss2_single_sample.wdl" as submit

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
    sample_id = assay_json['sample_id']

    fastq_1_name = assay_json['seq']['lanes'][0]['r1']
    fastq_2_name = assay_json['seq']['lanes'][0]['r2']
    fastq_1_url = name_to_meta[fastq_1_name]['url']
    fastq_2_url = name_to_meta[fastq_2_name]['url']

    print('Creating input map')
    with open('inputs.tsv', 'w') as f:
        f.write('fastq_1\tfastq_2\tsample_id\n')
        f.write('{0}\t{1}\t{2}\n'.format(fastq_1_url, fastq_2_url, sample_id))
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

workflow WrapperSs2RsemSingleSample {
  String bundle_uuid
  String bundle_version

  File gtf
  File ref_fasta
  File rrna_interval
  File ref_flat
  String star_genome
  String rsem_genome

  call GetInputs as prep {
    input:
      bundle_uuid = bundle_uuid,
      bundle_version = bundle_version
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

  call submit.SubmitSs2RsemSingleSample {
    input:
      bam_file = analysis.bam_file,
      bam_trans = analysis.bam_trans,
      rna_metrics = analysis.rna_metrics,
      aln_metrics = analysis.aln_metrics,
      rsem_gene_results = analysis.rsem_gene_results,
      rsem_isoform_results = analysis.rsem_isoform_results,
      rsem_gene_count = analysis.rsem_gene_count,
      gene_unique_counts = analysis.gene_unique_counts,
      exon_unique_counts = analysis.exon_unique_counts,
      transcript_unique_counts = analysis.transcript_unique_counts
  }
}
