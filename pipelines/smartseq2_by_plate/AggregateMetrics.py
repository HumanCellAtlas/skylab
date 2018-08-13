from crimson import picard
import loompy
import argparse
import os
import pandas as pd


def AggregatePicardMetricsRow(filenames, output_name):
    """
    pipeline output picard QC metrics at sinle cell/sample level.
    This functin is called to merge/aggregate QC metrics by metrics type and
    then merge multiple QC measurement into single matrix file.
    In this file, column is sample/cell and row is QC metrics
    :param files: metric files from pipeline outputs
    :param met_name: metrics name with workflow name and subworkflow name
    as prefix.such as 'run_pipelines.RunStarPipeline.alignment_summary_metrics'
    """
    # initial output
    mets = {}
    for file_name in filenames:
        cell_id = os.path.basename(file_name).split('_qc')[0]
        parsed = picard.parse(file_name)
        class_name = parsed['metrics']['class'].split('.')[2]
        # Aignment metrics return multiple lines,
        # but only output PAIRED-READS/third line
        if class_name == "AlignmentSummaryMetrics":
            # only parse out pair reads
            met = parsed['metrics']['contents'][2]
        # sometimes(very rare), insertion metrics also return multiple lines
        # results to include TANDEM repeats. but we only output the first line.
        elif class_name == "InsertSizeMetrics":
            # if the elemnet counts is less than 21,
            # it means insertion metrics returns multiple line results.
            if len(parsed['metrics']['contents']) < 21:
                met = parsed['metrics']['contents'][0]
            else:
                met = parsed['metrics']['contents']
        else:
            # other metrics(so far) only return one line results.
            met = parsed['metrics']['contents']
        mets[cell_id] = {
                k: met[k] for k in met if k not in
                ['SAMPLE', 'LIBRARY', 'READ_GROUP']
                }
    df = pd.DataFrame.from_dict(mets, orient='columns')
    df.insert(0, 'Class', class_name)
    df.to_csv(output_name + '.csv')


def AggregatePicardMetricsTable(filenames, output_name):
    df = pd.DataFrame()
    for file_name in filenames:
        cell_id = os.path.basename(file_name).split('_qc')[0]
        parsed = picard.parse(file_name)
        dat = pd.DataFrame.from_dict(parsed['metrics']['contents'])
        dat.insert(0, 'Sample', cell_id)
        df = df.append(dat, ignore_index=True)
    df.to_csv(output_name + '.csv')


def AggregateQCMetrics(filenames, output_name):
    df = pd.DataFrame()
    for file_name in filenames:
        dat = pd.read_csv(file_name)
        df = df.append(dat, ignore_index=True)
    df.rename(columns={'Unnamed: 0': 'Metrics'}, inplace=True)
    df.to_csv(output_name + '.csv')


def ParseHISATStats(filenames, output_name):
    mets = {}
    for file_name in filenames:
        if '_qc' in file_name:
            cell_id = os.path.basename(file_name).split('_qc')[0]
            tag = "HISAT2G"
        elif '_rsem' in file_name:
            cell_id = os.path.basename(file_name).split('_rsem')[0]
            tag = "HISAT2T"
        with open(file_name) as f:
            dat = f.readlines()
            d = [x.strip().split(':') for x in dat]
            d.pop(0)
            mets[cell_id] = {x[0]: x[1].strip().split(' ')[0] for x in d}
    df = pd.DataFrame.from_dict(mets, orient='columns')
    df.insert(0, "Class", tag)
    df.to_csv(output_name + '.csv')


def ParseRSENStats(filenames, output_name):
    mets = {}
    for file_name in filenames:
        cell_id = os.path.basename(file_name).split('_rsem')[0]
        i = 0
        with open(file_name) as f:
            while i < 3:
                if i == 0:
                    [N0, N1, N2, N_tot] = f.readline().strip().split(" ")
                elif i == 1:
                    [nUnique, nMulti, nUncertain] = \
                        f.readline().strip().split(" ")
                elif i == 2:
                    [nHits, read_type] = f.readline().strip().split(" ")
                i = i+1
        mets[cell_id] = {
                "unalignable reads": N0,
                "alignable reads": N1,
                "filtered reads": N2,
                "total reads": N_tot,
                "unique aligned": nUnique,
                "multiple mapped": nMulti,
                "total alignments": nHits,
                "strand": read_type,
                "uncertain reads": nUncertain
        }
    df = pd.DataFrame.from_dict(mets, orient='columns')
    df.insert(0, "Class", "RSEM")
    df.to_csv(output_name + '.csv')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--file_names",
        dest="file_names",
        nargs='+',
        required=True,
        help="a list of files to be parsed out."
    )
    parser.add_argument(
        "-o",
        "--output_name",
        dest="output_name",
        required=True,
        help="The output file name")
    parser.add_argument(
        "-t",
        "--metrics-type",
        dest="mtype",
        required=True,
        help="a list of string to represent metrics types,such Picard, PicardTable, HISAT2,RSEM, Core")
    args = parser.parse_args()
    if args.mtype == "Picard":
        AggregatePicardMetricsRow(args.file_names, args.output_name)
    elif args.mtype == "PicardTable":
        AggregatePicardMetricsTable(args.file_names, args.output_name)
    elif args.mtype == "Core":
        AggregateQCMetrics(args.file_names, args.output_name)
    elif args.mtype == "HISAT2":
        ParseHISATStats(args.file_names, args.output_name)
    elif args.mtype == "RSEM":
        ParseRSENStats(args.file_names, args.output_name)
    else:
        return 0


if __name__ == "__main__":
    main()
