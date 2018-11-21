import argparse
import csv
import re
import gzip
from scipy import sparse
import numpy as np
import zarr
from zarr import Blosc



ZARR_GROUP = {
    'expression_matrix': ["expression", "cell_id", "gene_id",
                          "gene_metadata_string", "gene_metadata_numeric", 
                           "cell_metadata_string", "cell_metadata_numeric"
                         ]
}

COMPRESSOR = Blosc(cname='lz4', clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)

# the number of rows in a chunk for expression counts
CHUNK_ROW_SIZE = 1000

def init_zarr(sample_id, path, file_format):
    """Initializes the zarr output
    Args:
        sample_id (str): sample or cell id
        path (str): path to the zarr output
        file_format (str): zarr file format [ DirectoryStore, ZipStore]
    Returns:
        zarr.hierarchy.Group
    """

    store = None
    if file_format == "DirectoryStore":
        store = zarr.DirectoryStore(path)

    if file_format == "ZipStore":
        store = zarr.ZipStore(path, mode='w')

    # create the root group
    root = zarr.group(store, overwrite=True)

    # add some readme for the user
    root.attrs['README'] = "The schema adopted in this zarr store may undergo changes in the the future"
    root.attrs['sample_id'] = sample_id

    # now iterate through list of expected groups and create them
    for dataset in ZARR_GROUP:
        root.create_group(dataset, overwrite=True)

    return root


def add_gene_metrics(data_group, input_path):
    """Converts the gene metrics from the Optimus pipeline to zarr file
    Args:
        data_group (zarr.hierarchy.Group): datagroup object for the zarr
        input_path (str): file containing gene metrics name and values
    """

    # read the gene metrics names and values
    if input_path.endswith(".gz"):
        with gzip.open(input_path, 'rt') as f:
            gene_metrics = [row for row in csv.reader(f)]
    else:
        with open(input_path, 'r') as f:
            gene_metrics = [row for row in csv.reader(f)]

    # metric names we use [1:] to remove the empty string
    if gene_metrics[0]:
        data_group.create_dataset(
            "gene_metadata_numeric_name",
            shape=(len(gene_metrics[0][1:]),),
            compressor=COMPRESSOR,
            dtype="<U40",
            chunks=(len(gene_metrics[0][1:]),),
            data=[list(gene_metrics[0][1:])])


    # Gene metric values, the row and column sizes
    ncols = len(gene_metrics[0][1:])
    nrows = len(gene_metrics[1:])

    # ignore the first line with the metric names in text
    gene_metric_values = []
    for row in gene_metrics[1:]:
        row_values = []
        for value_string in row[1:]:
            # some of the standard deviation values do not exist for one reads matches
            try:
                value = np.float32(value_string)
            except ValueError:
                value = np.nan
            row_values.append(value)
        gene_metric_values.append(row_values)

    # now insert the dataset that has the numeric values for the qc metrics for the genes
    data_group.create_dataset(
        "gene_metadata_numeric",
        shape=(nrows, ncols),
        compressor=COMPRESSOR,
        dtype=np.float32,
        chunks=(nrows, ncols),
        data=gene_metric_values)


def add_cell_metrics(data_group, input_path):
    """Converts cell metrics from the Optimus pipeline to zarr file
    Args:
        data_group (zarr.hierarchy.Group): datagroup object for the zarr
        input_path (str): file containing gene metrics name and values
    """

    # read the gene metrics names and values
    if input_path.endswith(".gz"):
        with gzip.open(input_path, 'rt') as f:
            cell_metrics = [row for row in csv.reader(f)]
    else:
        with open(input_path, 'r') as f:
            cell_metrics = [row for row in csv.reader(f)]

    # metric names for the cells
    data_group.create_dataset(
        "cell_metadata_numeric_name",
        shape=(len(cell_metrics[0][1:]), ),
        compressor=COMPRESSOR,
        dtype="<U40",
        chunks=(len(cell_metrics[0][1:]), ),
        data=[list(cell_metrics[0][1:])])

    # ignore the first line with the cell metric names in text
    cell_metric_values = []
    for row in cell_metrics[1:]:
        row_values = []
        for value_string in row[1:]:
            # some of the standard deviation values do not exist for one reads matches
            try:
                value = np.float32(value_string)
            except ValueError:
                value = np.nan
            row_values.append(value)
        cell_metric_values.append(row_values)

    # Gene metric values, the row and column sizes
    ncols = len(cell_metrics[0][1:])
    nrows = len(cell_metrics[1:])

    # now insert the dataset that has the numeric values for the cell qc metrics
    data_group.create_dataset(
        "cell_metadata_numeric",
        shape=(nrows, ncols),
        compressor=COMPRESSOR,
        dtype=np.float32,
        chunks=(nrows, ncols),
        data=cell_metric_values)


def add_expression_counts(data_group, args):
    """Converts  the count matrix from the Optimus pipeline to zarr file
    Args:
        data_group (zarr.hierarchy.Group): datagroup object for the zarr
        args (argparse.Namespace): input arguments for the run
    """
    # read the cell ids and adds into the cell_barcodes dataset
    barcodes = np.load(args.cell_ids)

    # note that if we do not specify the exact dimension, i.e., (len(barcodes), )
    # instead of (1, len(barcodes)) then  there is a memory bloat
    # while adding the cell_id or gene_id
    print('cell id')
    if len(barcodes):
        data_group.create_dataset(
            "cell_id",
            shape=(len(barcodes),),
            compressor=COMPRESSOR,
            dtype='<U40',
            #chunks=(len(barcodes), ),
            chunks=(10000,100),
            data=[list(barcodes)])

    # read the gene ids  and adds into the gene_ids dataset
    gene_ids = np.load(args.gene_ids)

    print('gene id')
    if len(gene_ids):
        data_group.create_dataset(
            "gene_id",
            shape=(len(gene_ids),),
            compressor=COMPRESSOR,
            dtype='<U40',
            #chunks=(len(gene_ids), ),
            chunks=(10000, 100),
            data=[list(gene_ids)])

    # read .npz file expression counts and add it to the expression_counts dataset
    exp_counts = np.load(args.count_matrix)
    # now convert it back to a csr_matrix object
    csr_exp_counts = sparse.csr_matrix((exp_counts['data'],
                                        exp_counts['indices'],
                                        exp_counts['indptr']),
                                       shape=exp_counts['shape'])

    print("exp", CHUNK_ROW_SIZE, csr_exp_counts.shape[1])
    # now create a dataset of zeros with the same dimensions as the expression count matrix
    exp_counts_group = data_group.zeros('expression',
                                        compressor=COMPRESSOR,
                                        shape=csr_exp_counts.shape,
                                        chunks=(CHUNK_ROW_SIZE, csr_exp_counts.shape[1]),
                                        dtype=np.float32)

    # load a chunks from the expression count matrix and update the corresponding chunk in expression matrix in zar
    for i in range(0, csr_exp_counts.shape[0], CHUNK_ROW_SIZE):
        # check if it is possible to make a full chunk of data
        if i + CHUNK_ROW_SIZE <= csr_exp_counts.shape[0]:
            exp_counts_group[i:i + CHUNK_ROW_SIZE, ] = csr_exp_counts[i:i + CHUNK_ROW_SIZE, ].toarray()
        else:
            # not enough data remaining, so make the final set of data less than a chunk
            j = csr_exp_counts.shape[0] - i
            exp_counts_group[i:i + j, ] = csr_exp_counts[i:i + j, ].toarray()


def create_zarr_files(args):
    """This function creates the zarr file or folder structure in output_zarr_path in format file_format,
        with sample_id from the input folder analysis_output_path
    Args:
        args (argparse.Namespace): input arguments for the run
    """
    # initiate the zarr file
    root_group = init_zarr(args.sample_id, args.output_zarr_path, args.zarr_format)

    # add the the gene metrics
    add_gene_metrics(root_group['expression_matrix'], args.gene_metrics)

    # add the the gene metrics
    add_cell_metrics(root_group['expression_matrix'], args.cell_metrics)

    # add the expression count matrix data
    print('hey')
    add_expression_counts(root_group['expression_matrix'], args)


def main():
    description = """This script converts the some of the Optimus outputs in to
                   zarr format (https://zarr.readthedocs.io/en/stable/) relevant output. 
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--cell_metrics',
                        dest="cell_metrics",
                        required=True,
                        help='a .csv file path for the cell metrics, an output of the MergeCellMetrics task')

    parser.add_argument('--gene_metrics',
                        dest="gene_metrics",
                        required=True,
                        help='a .csv file path for the gene metrics, an output of the MergeGeneMetrics task')

    parser.add_argument('--cell_id',
                        dest="cell_ids",
                        required=True,
                        help='a .npy file path for the cell barcodes, an output of the MergeCountFiles task')

    parser.add_argument('--gene_id',
                        dest="gene_ids",
                        required=True,
                        help='a .npy file path for the gene ids, an output of the MergeCountFiles task')

    parser.add_argument('--count_matrix',
                        dest="count_matrix",
                        required=True,
                        help='a .npz file path for the count matrix, an output of the MergeCountFiles task')

    parser.add_argument('--sample_id', dest="sample_id",
                        default="Unknown sample",
                        help='the sample name in the bundle')

    parser.add_argument('--output_path_for_zarr',
                        dest="output_zarr_path",
                        required=True,
                        help='path to .zarr file is to be created')

    parser.add_argument('--format', dest="zarr_format",
                        default="DirectoryStore",
                        choices=["DirectoryStore", "ZipStore"],
                        help='format of the zarr file choices: [DirectoryStore, ZipStore] default: DirectoryStore')
    args = parser.parse_args()

    create_zarr_files(args)


if __name__ == '__main__':
    main()
