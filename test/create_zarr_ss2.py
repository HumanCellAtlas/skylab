import argparse
import csv
import sys
import re, logging
import numpy
import zarr
from zarr import Blosc

""" key,  dataset name,   suffix of file,  metadata description  this is a list of groups to add to the zarr 
    file, which can be  further extended
"""

ZARR_GROUP = {
    'expression_matrix': ["qc_values", "qc_metric", "expression", "cell_id", "gene_id"]
}

COMPRESSOR = Blosc(cname='lz4', clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)


def init_zarr(sample_id, path, fileformat):
    """Initializes the zarr output
    Args:
        sample_id (str): sample or cell id
        path (str): path to the zarr output
        fileformat (str): zarr file format [ DirectoryStore, ZipStore]
    Returns:
        zarr.hierarchy.Group
    """

    store = None
    if fileformat == "DirectoryStore":
        store = zarr.DirectoryStore(path)

    if fileformat == "ZipStore":
        store = zarr.ZipStore(path, mode='w')

    # create the root group
    root = zarr.group(store, overwrite=True)

    # add some readme for the user
    root.attrs['README'] = "The schema adopted in this zarr store may undergo changes in future"
    root.attrs['sample_id'] = sample_id

    # now iterate through list of expected groups and create them
    for data in ZARR_GROUP.keys():
        datagroup = root.create_group(data, overwrite=True)

    return root


def read_and_convert_QC(datagroup, input_files_string):
    """Converts the QC of Smart Seq2 gene file pipeline outputs to zarr file
    Args:
        datagroup (zarr.hierarchy.Group): zarr group to add the datasets
        input_files_string (str): a string, comma separeted concatenated file names from the GroupQCOutputs task
                                    from the SmartSeq2SingleSample workflow
    """
    # split the files string by , 
    # create an array of the file names containing the suffix _QCs.csv
    relevant_qc_files = [x.strip() for x in input_files_string.split(',') if x.endswith('_QCs.csv')]

    # we expect only one, write an warning 
    if len(relevant_qc_files) != 1:
        logging.warning(sys.stderr.write("WARNING: Multiple files ending with _QCs.csv. Expected only one"))

    # read the QC values
    QC_values = []
    with open(relevant_qc_files[0], 'r') as f:
        QC_values = [row for row in csv.reader(f)]


    # Metrics
    # first value is the row name so we remove it"""
    datagroup.create_dataset(
        "qc_metric",
        shape=(1, len(QC_values[0][1:])),
        compressor=COMPRESSOR,
        dtype='<U40',
        chunks=(1, len(QC_values[0][1:])),
        data=[QC_values[0][1:]]
    )

    # Values, which are the values of the metrics
    datagroup.create_dataset(
        "qc_values",
        shape=(1, len(QC_values[2][1:])),
        compressor=COMPRESSOR,
        dtype="<U40",
        chunks=(1, len(QC_values[2][1:])),
        data=[QC_values[2][1:]]
    )

    # Cell IDs
    datagroup.create_dataset(
        "cell_id",
        shape=(1, len(QC_values[2][0:1])),
        compressor=COMPRESSOR,
        dtype="<U40",
        chunks=(1, len(QC_values[2][0:1])),
        data=[QC_values[2][0:1]]
    )


def read_and_convert_expression(datagroup, input_path):
    """Converts the Smart Seq2 gene file pipeline outputs to zarr file
    Args:
        datagroup (zarr.hierarchy.Group): databroup object for the zarr
        input_path (str): file where the SS2 pipeline expression counts are
    """
    reader = csv.DictReader(open(input_path), delimiter="\t")

    expression_values = {}
    for row in reader:
        expression_values[row["gene_id"]] = float(row["TPM"])

    # TPM
    datagroup.create_dataset(
        "expression",
        shape=(1, len(expression_values)),
        compressor=COMPRESSOR,
        dtype=numpy.float32,
        chunks=(1, len(expression_values)),
        data=[list(expression_values.values())]
    )

    # Gene IDs
    datagroup.create_dataset(
        "gene_id",
        shape=(1, len(expression_values)),
        compressor=COMPRESSOR,
        dtype="<U40",
        chunks=(1, len(expression_values)),
        data=[list(expression_values.keys())]
    )


def create_zarr_files(sample_id, qc_analysis_output_files_string, rsem_genes_results_file, output_zarr_path,
                      fileformat = "DirectoryStore"):
    """This function creates the zarr file or folder structure in output_zarr_path in format fileformat,
        with sample_id from the input folder analysis_output_path
    Args:
        sample_id (str): sample or cell id
        qc_analysis_output_files_string (str): a string with the file names in the QCGroup of SS2 pipeline output, separated by commas
        rsem_genes_results_file (str): the file for the expression count
        output_zarr_path (str): location of the output zarr
        fileformat (str): zarr file format  [DirectoryStore or ZipStore ] Default: DirectoryStore
    """
    # initiate the zarr file
    rootgroup = init_zarr(sample_id, output_zarr_path, fileformat)

    # add the the gene/TMP counts
    read_and_convert_expression(rootgroup['expression_matrix'], rsem_genes_results_file)

    # add the the QC metrics
    read_and_convert_QC(rootgroup['expression_matrix'], qc_analysis_output_files_string)


def main():
    description = """This script converts the some of the SmartSeq2 pipeline outputs in to
                   zarr format (https://zarr.readthedocs.io/en/stable/) relevant output. 
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--qc_analysis_output_files_string',
                        dest="qc_analysis_output_files_string",
                        help='a string, list of files from the GroupQCOutputs task of SS2 Sin Samp workflow')

    parser.add_argument('--rsem_genes_results',
                        dest="rsem_genes_results",
                        help='path to the folder containing the files to be added to the zarr')

    parser.add_argument('--output_path_for_zarr',
                        dest="output_zarr_path",
                        help='path where the zarr file is to be created')

    parser.add_argument('--sample_id', dest="sample_id",
                        default="Unknown sample",
                        help='the sample name in the bundle')

    parser.add_argument('--format', dest="format",
                        default="DirectoryStore",
                        help='format of the zarr file choices: [DirectoryStore, ZipStore] default: DirectoryStore')
    args = parser.parse_args()

    create_zarr_files(args.sample_id, args.qc_analysis_output_files_string, args.rsem_genes_results,
                      args.output_zarr_path, args.format)


if __name__ == '__main__':
    main()
