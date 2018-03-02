import argparse
import pandas as pd
import numpy as np
from google.cloud import storage
import json
from os.path import basename
import sys
import requests


def retrieve_workflow_outputs(uuid, run_name):
    # load cromwell credential
    logins = json.load(open('/usr/secrets/broad-dsde-mint-dev-cromwell.json'))
    metadata_url = "https://cromwell.mint-dev.broadinstitute.org/api/workflows/v1/" + uuid + "/metadata?expandSubWorkflows=false"
    r = requests.get(
        metadata_url,
        auth=(logins['cromwell_username'], logins['cromwell_password']))
    data = r.json()
    # load output files
    files = data['outputs'][run_name]
    return (files)


def counts_to_tpm(counts, lengths):
    rpk = [float(s) / float(l) for s, l in zip(counts, lengths)]
    RPKM = sum(rpk) / 1e6
    TPM = [r / RPKM for r in rpk]
    return (TPM)


def merge_featurecount_outputs(files, value_name):
    """
    pipeline produce gene quantification by using either FeatureCounts at single sample level 
    This function is called to merge or aggregate all the single cell/sample level results into
    a single matrix file. FeatureCount only produce counts data so if value_name is tpm, first
    we need to convert counts into tpm
    param files: the output files(quantification) of pipeline run.
    param value_name: aggregate value type, can be est_count or tpm.
    """
    client = storage.Client()
    bucket = client.get_bucket('broad-dsde-mint-dev-cromwell-execution')
    # initial dataframe
    merged = pd.DataFrame()
    # loop through input files
    for kk in range(0, len(files)):
        fc1 = files[kk]
        fc1 = fc1.replace('gs://broad-dsde-mint-dev-cromwell-execution/', '')
        blob1 = bucket.get_blob(fc1)
        bname1 = basename(fc1)
        # sample_name is the prefix of input file name, such as SRR123456
        sample_name = bname1.split('.')[0]
        with open(bname1, 'wb') as file_obj:
            blob1.download_to_file(file_obj)
        # load file and only select gene_id, lenght and sample column
        dat = pd.read_csv(
            bname1,
            skiprows=2,
            sep='\t',
            usecols=['gene_id', 'length', sample_name],
            names=[
                'gene_id', 'contigs', 'start', 'end', 'strands', 'length',
                sample_name
            ])
        # featurecount only output counts so need to convert counts to tpm
        if value_name == 'tpm':
            sample = list(dat[sample_name].values)
            lengths = list(dat['length'].values)
            TPM = counts_to_tpm(sample, lengths)
            dat[sample_name] = TPM
        # while looping input files, merging parsed results into a single file as output
        # initial merged file by the first input file
        if kk == 0:
            merged = dat
        else:
            # only parse out gene_id and sample column from input files if it is not the first input file.
            cnt = dat
            # merge by gene_id
            merged = pd.merge(
                left=merged,
                right=cnt[['gene_id', sample_name]],
                left_on='gene_id',
                right_on='gene_id')
    merged = merged.round(3)
    return (merged)


def merge_rsem_outputs(files, value_name):
    """
    pipeline produce gene quantification by using either RSEM at single sample level 
    This function is called to merge or aggregate all the single cell/sample level results into
    a single matrix file
    param files: the output files(quantification) of pipeline run.
    param value_name: aggregate value type, can be est_count or tpm. 
    """
    # set up auth
    client = storage.Client()
    bucket = client.get_bucket('broad-dsde-mint-dev-cromwell-execution')
    # initial dataframe
    merged = pd.DataFrame()
    # loop through input files
    for kk in range(0, len(files)):
        fc1 = files[kk]
        fc1 = fc1.replace('gs://broad-dsde-mint-dev-cromwell-execution/', '')
        blob1 = bucket.get_blob(fc1)
        bname1 = basename(fc1)
        # sample_name is the prefix of input file name, such as SRR123456
        sample_name = bname1.split('.')[0]
        with open(bname1, 'wb') as file_obj:
            blob1.download_to_file(file_obj)
        # if it is first iteration, need to parse out header
        if kk == 0:
            merged = pd.read_csv(
                bname1,
                delimiter='\t',
                skiprows=1,
                sep='\t',
                names=[
                    'gene_id', 'trans_id', 'length', 'eff_length',
                    'est_counts', 'tpm', 'fpkm', 'pme_count', 'pme_sd',
                    'pme_tpm', 'pme_fpkm'
                ],
                usecols=['gene_id', 'length', value_name])
            merged.columns = ['gene_id', 'length', sample_name]
        else:
            cnt = pd.read_csv(
                bname1,
                delimiter='\t',
                skiprows=1,
                sep='\t',
                names=[
                    'gene_id', 'trans_id', 'length', 'eff_length',
                    'est_counts', 'tpm', 'fpkm', 'pme_count', 'pme_sd',
                    'pme_tpm', 'pme_fpkm'
                ],
                usecols=['gene_id', value_name])
            cnt.columns = ['gene_id', sample_name]
            merged = pd.merge(
                left=merged, right=cnt, left_on='gene_id', right_on='gene_id')
    # round matrix by 3 digits
    merged = merged.round(3)
    return (merged)


def run_merge_expression(uuid, run_name, output_name, value_name):
    """
    this function will call a parsing function to merge expression vectors 
    based on run_name, for example, if run_name is unq_gene_counts or mult_genes_counts,
    it indicate a FeatureCount outputs
    if rsem in run_name, it indicates  rsem outputs
    """
    # retrieve workflows outputs
    vfiles = retrieve_workflow_outputs(uuid, run_name)
    # check run_time and then decide which function to call
    if "unq_genes_counts" in run_name or "mult_genes_counts" in run_name:
        merged_matrix = merge_featurecount_outputs(vfiles, value_name)
    if 'rsem' in run_name:
        merged_matrix = merge_rsem_outputs(vfiles, value_name)
    merged_matrix.to_csv(output_name, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-u",
        "--uuid",
        dest="uuid",
        required=True,
        help="The uuid of workflow, top level")
    parser.add_argument(
        "-rn",
        "--run_name",
        dest="run_name",
        required=True,
        help=
        "The output file name from task or workflow,ex RunSTARPipeline.rsem_gene_results"
    )
    parser.add_argument(
        "-o",
        "--output_name",
        dest="output_name",
        required=True,
        help="The output folder folder name.")
    parser.add_argument(
        "-t",
        "--value_name",
        dest="value_name",
        required=True,
        help="value or column to parse, ex, tpm, est_counts")
    args = parser.parse_args()
    run_merge_expression(args.uuid, args.run_name, args.output_name,
                         args.value_name)


if __name__ == "__main__":
    main()
