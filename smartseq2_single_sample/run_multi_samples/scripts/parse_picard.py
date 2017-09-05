from crimson import picard
import json
import pandas as pd
import argparse
import sys
from collections import OrderedDict
import os

def parse_picard_to_table_by_sample(input_metrics_list,output_filename):
    mets=OrderedDict()
    if os.stat(input_metrics_list).st_size==0:
        print("Input is empty.")
        sys.exit()
    with open(input_metrics_list) as f:
        for line in f:
            met_fn = line.strip('\n')
            base = os.path.basename(met_fn)
            [sra_id,ftype] = os.path.splitext(base)
            parsed = picard.parse(met_fn)
            class_name = parsed['metrics']['class']
            if class_name == "picard.analysis.AlignmentSummaryMetrics":
                met = parsed['metrics']['contents'][2]
            else:
                met = parsed['metrics']['contents']
            mets.update(met)

    tab = pd.DataFrame.from_records([mets])
    tab.index=[sra_id]
    tab.to_csv(output_filename+'.metrics')

def parse_picard_to_table_by_type(input_metrics_list,output_filename):
    mets={}
    with open(input_metrics_list) as f:
        for line in f:
            met_fn = line.strip('\n')
            base = os.path.basename(met_fn)
            [sra_id,ftype] = os.path.splitext(base)
            parsed = picard.parse(met_fn)
            class_name = parsed['metrics']['class']
            if class_name == "picard.analysis.AlignmentSummaryMetrics":
                met = parsed['metrics']['contents'][2]
            else:
                met = parsed['metrics']['contents']
            mets[sra_id] = met
    f.close()
    tab = pd.DataFrame.from_dict(mets)
    tab.to_csv(output_filename+'.metrics')

def merge_sample_to_table(input_metrics_list,output_filename):
    mets = pd.DataFrame()
    with open(input_metrics_list) as f:
        for line in f:
            met_fn = line.strip('\n')
            tab=pd.read_csv(met_fn)
            mets=mets.append(tab)
    f.close()
    mets.to_csv(output_filename+'_metrics.csv',index=None)

def main():
       
    parser = argparse.ArgumentParser()
    ## list of input file, could be picard metrics, matrix table.
    parser.add_argument("-I","--input-file",dest="input_metrics_file",required=True,help="list of metrics file to parse")
    ## action to take place: if type is 'sample',which means inputs are picard metrics of one sample, this script will aggregate all input metrics into one table. if type is  'metrics', which means the inputs are metrics of one picard collection, such as rna, aln, this script will aggregate all metrics into one table. if option is 'table', which means inputs are table matrix, possiblely be the outputs of 'sample', this scripts will append all input tables into one table. 
    parser.add_argument("-T","--type",dest="aggregate_type",required=True,help="aggregated by: sample or metrics or table.'sample' option merge several metrics files from one sample to single table.'metrics' option aggregates metrics files into a table by metrics type, such as rna_metrics, aln_metrics. 'table' option simplely append each input metrics table into a single table file.")
    parser.add_argument('-O',"--output",dest="output_filename",required=True,help="output file name or prefix name")
    args = parser.parse_args()
    if args.aggregate_type == "sample":
        parse_picard_to_table_by_sample(args.input_metrics_file,args.output_filename)
    elif args.aggregate_type == "metrics":
        parse_picard_to_table_by_type(args.input_metrics_file,args.output_filename)
    elif args.aggregate_type == "table":
        merge_sample_to_table(args.input_metrics_file,args.output_filename)
    else:
        print('Wrong input values for -T')
        sys.exit()

if __name__ == "__main__":
    main()
   
