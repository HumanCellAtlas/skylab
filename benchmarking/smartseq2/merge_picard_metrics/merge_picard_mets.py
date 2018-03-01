from crimson import picard
import pandas as pd
import numpy as np
from google.cloud import storage
import json
from os.path import basename
import sys
import requests
import argparse


def retrieve_workflow_outputs(uuid, run_name):
    # load cromwell credential
    logins = json.load(open('/usr/secrets/broad-dsde-mint-dev-cromwell.json'))
    # meta_url
    metadata_url = "https://cromwell.mint-dev.broadinstitute.org/api/workflows/v1/" + uuid + "/metadata?expandSubWorkflows=false"
    r = requests.get(
        metadata_url,
        auth=(logins['cromwell_username'], logins['cromwell_password']))
    data = r.json()
    # load output files
    files = data['outputs'][run_name]
    return (files)


def merge_picard_metrics(files, metric_name):
    """
    piepline output picard QC metrics at sinle cell/sample level.
    This functin is called to merge/aggregate QC metrics by metrics type and then merge multiple QC measurement 
    into single matrix file. In this file, column is sample/cell and row is QC metrics
    param files: metric files from pipeline outputs
    param met_name: metrics name with workflow name and subworkflow name as prefix. such as 'run_pipelines.RunStarPipeline.alignment_summary_metrics'
    """
    # set up auth
    client = storage.Client()
    bucket = client.get_bucket('broad-dsde-mint-dev-cromwell-execution')
    # load cromwell credential
    logins = json.load(open('/usr/secrets/broad-dsde-mint-dev-cromwell.json'))
    # initial output
    mets = {}
    for kk in range(0, len(files)):
        fc1 = files[kk]
        fc1 = fc1.replace('gs://broad-dsde-mint-dev-cromwell-execution/', '')
        blob1 = bucket.get_blob(fc1)
        bname1 = basename(fc1)
        # sample name is prefix of file name
        sample_name = bname1.split('.')[0]
        with open(bname1, 'wb') as file_obj:
            blob1.download_to_file(file_obj)
        # use picard package parse out picard output, a json file is returned
        parsed = picard.parse(bname1)
        class_name = parsed['metrics']['class']
        # Aignment metrics return multiple lines, but only output PAIRED-READS/third line
        if class_name == "picard.analysis.AlignmentSummaryMetrics":
            ## only parse out pair reads
            met = parsed['metrics']['contents'][2]
        # sometimes(very rare), insertion metrics also return multiple lines results to include TANDEM repeats. but we only output the first line.
        elif class_name == "picard.analysis.InsertSizeMetrics":
            # if the elemnet counts is less than 21, it means insertion metrics returns multiple line results.
            if len(parsed['metrics']['contents']) < 21:
                met = parsed['metrics']['contents'][0]
            else:
                met = parsed['metrics']['contents']
        else:
            # other metrics(so far) only return one line results.
            met = parsed['metrics']['contents']
        mets[sample_name] = met
    tab = pd.DataFrame.from_dict(mets)
    return (tab)


def run_merge_metrics(uuid, metric_name, output_name):
    """
    call functions to nerge metrics and output in one file
    """
    metfiles = retrieve_workflow_outputs(uuid, metric_name)
    metrics_matrix = merge_picard_metrics(metfiles, metric_name)
    metrics_matrix.to_csv(output_name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-u",
        "--uuid",
        dest="uuid",
        required=True,
        help="The uuid of workflow, top level")
    parser.add_argument(
        "-m",
        "--metrics_name",
        dest="met_name",
        required=True,
        help="The Picard metrics class name")
    parser.add_argument(
        "-o",
        "--output_name",
        dest="output_name",
        required=True,
        help="The output file name")
    args = parser.parse_args()
    run_merge_metrics(args.uuid, args.met_name, args.output_name)


if __name__ == "__main__":
    main()
