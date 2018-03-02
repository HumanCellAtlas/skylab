import argparse
import pandas as pd
import numpy as np
from google.cloud import storage
import json
from os.path import basename
import sys
import requests


def retrieve_workflow_outputs(cromwell_uuid, file_name):
    # load cromwell credential
    logins = json.load(open('/usr/secrets/broad-dsde-mint-dev-cromwell.json'))
    metadata_url = "https://cromwell.mint-dev.broadinstitute.org/api/workflows/v1/" + cromwell_uuid + "/metadata?expandSubWorkflows=false"
    r = requests.get(
        metadata_url,
        auth=(logins['cromwell_username'], logins['cromwell_password']))
    data = r.json()
    # load output files
    files = data['outputs'][file_name]
    return files


def counts_to_tpm(counts, lengths):
    rpk = [float(s) / float(l) for s, l in zip(counts, lengths)]
    rpkm = sum(rpk) / 1e6
    tpm = [r / rpkm for r in rpk]
    return tpm


def merge_featurecount_outputs(files, value_name):
    """
    pipeline produce gene quantification by using either FeatureCounts at single sample level 
    This function is called to merge or aggregate all the single cell/sample level results into
    a single matrix file. FeatureCount only produce counts data so if value_name is tpm, first
    we need to convert counts into tpm
    :param files: the output files(quantification) of pipeline run.
    :param value_name: aggregate value type, can be est_count or tpm.
    :return function returns a dataframe *merged*  which include matrix of counts/tpm, 
    row represent gene and column represents 
    sample/cell
    """
    client = storage.Client()
    bucket = client.get_bucket('broad-dsde-mint-dev-cromwell-execution')
    # initial dataframe
    merged = pd.DataFrame()
    # loop through input files
    for kk in range(0, len(files)):
        fc = files[kk]
        fc = fc.replace('gs://broad-dsde-mint-dev-cromwell-execution/', '')
        blob = bucket.get_blob(fc)
        exp_name = basename(fc)
        # sample_name is the prefix of input file name, such as SRR123456
        sample_name = exp_name.split('.')[0]
        with open(exp_name, 'wb') as file_obj:
            blob.download_to_file(file_obj)
        # load file and only select gene_id, length and sample column
        dat = pd.read_csv(
            exp_name,
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
            tpm = counts_to_tpm(sample, lengths)
            dat[sample_name] = tpm
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
    return merged


def merge_rsem_outputs(files, value_name):
    """
    pipeline produce gene quantification by using either RSEM at single sample level 
    This function is called to merge or aggregate all the single cell/sample level results into
    a single matrix file
    :param files: the output files(quantification) of pipeline run.
    :param value_name: aggregate value type, can be est_count or tpm. 
    :return function returns a dataframe *merged*  which include matrix of counts/tpm,
    row represent gene and column represents 
    """
    # set up auth
    client = storage.Client()
    bucket = client.get_bucket('broad-dsde-mint-dev-cromwell-execution')
    # initial dataframe
    merged = pd.DataFrame()
    # loop through input files
    for kk in range(0, len(files)):
        fc = files[kk]
        fc = fc.replace('gs://broad-dsde-mint-dev-cromwell-execution/', '')
        blob = bucket.get_blob(fc)
        exp_name = basename(fc)
        # sample_name is the prefix of input file name, such as SRR123456
        sample_name = exp_name.split('.')[0]
        with open(exp_name, 'wb') as file_obj:
            blob.download_to_file(file_obj)
        # laod data, only select gene id, length and a data column
        dat = pd.read_csv(
            exp_name,
            delimiter='\t',
            skiprows=1,
            sep='\t',
            names=[
                'gene_id', 'trans_id', 'length', 'eff_length',
                'est_counts', 'tpm', 'fpkm', 'pme_count', 'pme_sd',
                'pme_tpm', 'pme_fpkm'
                ],  
            usecols=['gene_id', 'length', value_name])
        # while looping input files, merging parsed results into a single file as output
        # initial merged file by the first input file
        dat.columns = ['gene_id','length',sample_name]
        if kk == 0:
            merged = dat
        else:
            cnt = dat
            merged = pd.merge(
                left=merged, right=cnt[['gene_id',sample_name]], left_on='gene_id', right_on='gene_id')
    # round matrix by 3 digits
    merged = merged.round(3)
    return merged


def run_merge_expression(cromwell_uuid, file_name, output_name, value_name):
    """
    this function will call a parsing function to merge expression vectors 
    :param cromwell_uuid cromwell workflow uuid
    :param file_name can be rsem or mult_gene_results for featurecounts
    :param value_name can be est_counts or tpm
    :param output_name output csv file name
    """
    # retrieve workflows outputs
    vfiles = retrieve_workflow_outputs(cromwell_uuid, file_name)
    # check run_time and then decide which function to call
    if "unq_genes_counts" in file_name or "mult_genes_counts" in file_name:
        merged_matrix = merge_featurecount_outputs(vfiles, value_name)
    if 'rsem' in file_name:
        merged_matrix = merge_rsem_outputs(vfiles, value_name)
    merged_matrix.to_csv(output_name, index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-u",
        "--uuid",
        dest="uuid",
        required=True,
        help="The uuid of main workflow")
    parser.add_argument(
        "-fn",
        "--file_name",
        dest="file_name",
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
    run_merge_expression(args.uuid, args.file_name, args.output_name,
                         args.value_name)


if __name__ == "__main__":
    main()
