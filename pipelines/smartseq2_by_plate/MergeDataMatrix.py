import pandas as pd
from os.path import basename
import argparse

def MergeRsemQuantification(files, col_name, output):
    """
    pipeline produce gene quantification by using either RSEM at single sample level 
    This function is called to merge or aggregate all the single cell/sample level results into
    a single matrix file
    :param files: the list of rsem outputs.
    :param col_name: aggregate value type, can be est_counts or tpm. 
    :param utput: output file name, *,csv or *.tsv
    """
    
    # initial dataframe
    merged = pd.DataFrame()
    filenames = files.split(',')
    # loop through input files
    for kk in range(0, len(filenames)):
        file_name = filenames[kk]
        base_name = basename(file_name)
        # sample_name is the prefix of input file name, such as SRR123456
        sample_name = base_name.split('.')[0].split('_')[0]
        btype = base_name.split('.')[1]
        # load data, only select gene id, length and a data column
        if btype == "genes":
            dat = pd.read_csv(
                file_name,
                delimiter='\t',
                skiprows=1,
                sep='\t',
                names=[
                    'gene_id', 'trans_id', 'length', 'eff_length',
                    'est_counts', 'tpm', 'fpkm', 'pme_count', 'pme_sd',
                    'pme_tpm', 'pme_fpkm'
                    ],
                        usecols=['gene_id', 'length', col_name])
        # while looping input files, merging parsed results into a single file as output
        # initial merged file by the first input file
        else:
            dat = pd.read_csv(
                file_name,
                delimiter='\t',
                skiprows=1,
                sep='\t',
                names=[
                    'trans_id','gene_id','length','eff_length',
                    'est_counts','tpm','fpkm','IsoPct','pm_count',
                    'pme_std','pme_tpm','pme_fpkm','IsoPct_pme_tpm'
                    ],
                        usecols=['trans_id', 'length', col_name])
        # while looping input files, merging parsed results into a single file as output
        # initial merged file by the first input file
        dat.columns = ['tracking_id','length',sample_name]
        if kk == 0:
            merged = dat
        else:
            cnt = dat
            merged = pd.merge(
                left=merged, right=cnt[['tracking_id',sample_name]],
                left_on='tracking_id', right_on='tracking_id')
    # round matrix by 3 digits
    merged = merged.round(3)
    merged.to_csv(output, index=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--file_names",
        dest="file_names",
        required=True,
        help=
        "a list of files to be parsed out."
    )
    parser.add_argument(
        "-o",
        "--output_name",
        dest="output_name",
        required=True,
        help="The output file name")
    parser.add_argument(
        "-t",
        "--col_name",
        dest="col_name",
        required=True,
        help="name or column to parse, possible tpm, est_counts")
    args = parser.parse_args()
    MergeRsemQuantification(args.file_names, args.col_name,args.output_name)

if __name__ == "__main__":
    main()
