import json
import argparse
import pandas as pd

def AggregatePicardMetric(filenames, output_name):
    """
    piepline output picard QC metrics at sinle cell/sample level.
    This functin is called to merge/aggregate QC metrics by metrics type and then merge multiple QC measurement 
    into single matrix file. In this file, column is sample/cell and row is QC metrics
    :param files: metric files from pipeline outputs
    :param output_name: output file name
    """
    # initial output
    print(filenames)
    for kk in range(0, len(filenames)):
        d = json.load(open(filenames[kk]))
        if kk == 0 :
            mets=d.copy()
        else:
            mets.update(d)
    merged = pd.DataFrame.from_dict(mets)
    merged.to_csv(output_name+'.picard.batch.core.csv')
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--file_names",
        dest="file_names",
        required=True,
        nargs='+',
        help=
        "a list of files to be parsed out."
    )
    parser.add_argument(
        "-o",
        "--output_name",
        dest="output_name",
        required=True,
        help="The output file name")
    args = parser.parse_args()
    AggregatePicardMetric(args.file_names, args.output_name)

if __name__ == "__main__":
    main()
