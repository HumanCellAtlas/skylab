import argparse
import csv
import gzip
import numpy as np
import zarr
from scipy import sparse
from zarr import Blosc
import logging

ZARR_GROUP = {
    'expression_matrix': ["expression", "cell_id", "gene_id",
                          "gene_metadata_string", "gene_metadata_numeric",
                          "cell_metadata_string", "cell_metadata_numeric"
                          ]
}

COMPRESSOR = Blosc(cname='lz4', clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)

# the number of rows in a chunk for expression counts
CHUNK_ROW_SIZE = 10000
CHUNK_COL_SIZE = 10000
logging.basicConfig(level=logging.INFO)


def init_zarr(sample_id, path, file_format):
    """Initializes the zarr output.

    Args:
        sample_id (str): sample or cell id
        path (str): path to the zarr output
        file_format (str): zarr file format [DirectoryStore, ZipStore]

    Returns:
        root (zarr.hierarchy.Group): initialized zarr group
    """

    store = None
    if file_format == "DirectoryStore":
        store = zarr.DirectoryStore(path)

    if file_format == "ZipStore":
        store = zarr.ZipStore(path, mode='w')

    # create the root group
    root = zarr.group(store, overwrite=True)

    # add some readme for the user
    root.attrs['README'] = "The schema adopted in this zarr store may undergo changes in the future"
    root.attrs['sample_id'] = sample_id

    # now iterate through list of expected groups and create them
    for dataset in ZARR_GROUP:
        root.create_group(dataset, overwrite=True)

    return root


def add_gene_metrics(data_group, input_path, gene_ids, verbose=False):
    """Converts the gene metrics from the Optimus pipeline to zarr file

    Args:
        data_group (zarr.hierarchy.Group): datagroup object for the zarr
        input_path (str): file containing gene metrics name and values
        gene_ids (list): list of gene ids
        verbose (bool): whether to output verbose messages for debugging purposes
    """

    # read the gene metrics names and values
    if input_path.endswith(".gz"):
        with gzip.open(input_path, 'rt') as f:
            gene_metrics = [row for row in csv.reader(f)]
    else:
        with open(input_path, 'r') as f:
            gene_metrics = [row for row in csv.reader(f)]

    # metric names we use [1:] to remove the empty string
    if len(gene_metrics[0][1:]):
        data_group.create_dataset(
            "gene_metadata_numeric_name",
            shape=(len(gene_metrics[0][1:]),),
            compressor=COMPRESSOR,
            dtype="<U40",
            chunks=(len(gene_metrics[0][1:]),),
            data=list(gene_metrics[0][1:]))
    else:
        logging.info("Not adding \"gene_metadata_numeric_name\" to zarr output: must have at least one metric")

    if verbose:
        logging.info("# gene numeric metadata", len(gene_metrics[0][1:]))

    # Gene metric values, the row and column sizes
    gene_ids_location = {gene_id: index for index, gene_id in enumerate(gene_ids)}

    # ignore the first line with the metric names in text
    ncols = 0
    gene_id_to_metric_values = {}
    for row in gene_metrics:
        # only consider genes that are in the count matrix
        if not row[0] in gene_ids_location:
            continue

        row_values = []
        for value_string in row[1:]:
            # some of the standard deviation values do not exist for one reads matches
            try:
                value = np.float32(value_string)
            except ValueError:
                value = np.nan
            row_values.append(value)
        gene_id_to_metric_values[row[0]] = row_values

        # note that all of these lengths are assumed to be equal and this check is already done in the pipeline
        ncols = len(row_values)

    # now insert the metrics of the cells that are in count matrix, i.e.,  the global variable "cell_ids"

    gene_metric_values = []
    for gene_id in gene_ids:
        if gene_id in gene_id_to_metric_values:
            gene_metric_values.append(gene_id_to_metric_values[gene_id])
        else:
            # if no metrics for a cell present in the count matrix then fill them with np.nans
            gene_metric_values.append([np.nan] * ncols)

    nrows = len(gene_ids)
    if verbose:
        logging.info("# of genes: {}".format(nrows))
        logging.info("# of gene metadate metrics: {}".format(ncols))

    # now insert the dataset that has the numeric values for the qc metrics for the genes
    if nrows and ncols:
        data_group.create_dataset(
            "gene_metadata_numeric",
            shape=(nrows, ncols),
            compressor=COMPRESSOR,
            dtype=np.float32,
            chunks=(nrows, ncols),
            data=gene_metric_values)
    else:
        logging.info("Not adding \"gene_metadata_numeric\" to zarr output: either the #genes or # cell ids is 0")


def add_cell_metrics(data_group, input_path, cell_ids, verbose=False):
    """Converts cell metrics from the Optimus pipeline to zarr file

    Args:
        data_group (zarr.hierarchy.Group): datagroup object for the zarr
        input_path (str): file containing gene metrics name and values
        cell_ids (list): list of cell ids
        verbose (bool): whether to output verbose messages for debugging purposes
    """
    # read the gene metrics names and values
    if input_path.endswith(".gz"):
        with gzip.open(input_path, 'rt') as f:
            cell_metrics = [row for row in csv.reader(f)]
    else:
        with open(input_path, 'r') as f:
            cell_metrics = [row for row in csv.reader(f)]

    # metric names for the cells
    if len(cell_metrics[0][1:]):
        data_group.create_dataset(
            "cell_metadata_numeric_name",
            shape=(len(cell_metrics[0][1:]),),
            compressor=COMPRESSOR,
            dtype="<U40",
            chunks=(len(cell_metrics[0][1:]),),
            data=list(cell_metrics[0][1:]))
    else:
        logging.info("Not adding \"cell_metadata_numeric_name\" to zarr output: must have at least one metric")

    if verbose:
        logging.info("# of cell_metadata_numeric_names: {}".format(len(cell_metrics[0][1:])))

    cell_ids_location = {cell_id: index for index, cell_id in enumerate(cell_ids)}

    cell_id_to_metric_values = {}
    # ignore the first line with the cell metric names in text

    ncols = 0
    for row in cell_metrics[1:]:
        # only consider cell_id that are also in the count_metrics
        if not row[0] in cell_ids_location:
            continue
        row_values = []
        for value_string in row[1:]:
            # some of the standard deviation values do not exist for one reads matches
            try:
                value = np.float32(value_string)
            except ValueError:
                value = np.nan
            row_values.append(value)
        cell_id_to_metric_values[row[0]] = row_values

        # note that all of these lengths are assumed to be equal and this check is already done in the pipeline
        if ncols == 0:
            ncols = len(row_values)

    # now insert the metrics of the cells that are in count matrix, i.e.,  the global variable "cell_ids"
    cell_metric_values = []
    for cell_id in cell_ids:
        if cell_id in cell_id_to_metric_values:
            cell_metric_values.append(cell_id_to_metric_values[cell_id])
        else:
            # if no metrics for a cell present in the count matrix then fill them with np.nans
            cell_metric_values.append([np.nan] * ncols)

    # Gene metric values, the row size (i.e., number of cell_ids)
    nrows = len(cell_ids)
    if verbose:
        logging.info('# of cell ids: {}'.format(nrows))
        logging.info('# of numeric metrics: {}'.format(ncols))

    # now insert the dataset that has the numeric values for the cell qc metrics
    if nrows and ncols:
        data_group.create_dataset(
            "cell_metadata_numeric",
            shape=(nrows, ncols),
            compressor=COMPRESSOR,
            dtype=np.float32,
            chunks=(nrows, ncols),
            data=cell_metric_values)
    else:
        logging.info("Not adding \"cell_metadata_numeric\" to zarr output: either the #genes or # cell ids is 0")


def add_expression_counts(data_group, args):
    """Converts  the count matrix from the Optimus pipeline to zarr file

    Args:
        data_group (zarr.hierarchy.Group): datagroup object for the zarr
        args (argparse.Namespace): input arguments for the run

    Return:
        cell_ids: list of cell ids
        gene_ids: list of gene ids
    """
    # read the cell ids and adds into the cell_barcodes dataset
    cell_ids = np.load(args.cell_ids)

    # note that if we do not specify the exact dimension, i.e., (len(barcodes), )
    # instead of (1, len(barcodes)) then  there is a memory bloat
    # while adding the cell_id or gene_id
    if len(cell_ids):
        data_group.create_dataset(
            "cell_id",
            shape=(len(cell_ids),),
            compressor=COMPRESSOR,
            dtype='<U40',
            chunks=(10000,),
            data=list(cell_ids))
    else:
        logging.info("Not adding \"cell_id\" to zarr output: # cell ids is 0")

        # read the gene ids  and adds into the gene_ids dataset
    gene_ids = np.load(args.gene_ids)

    if len(gene_ids):
        data_group.create_dataset(
            "gene_id",
            shape=(len(gene_ids),),
            compressor=COMPRESSOR,
            dtype='<U40',
            chunks=(10000,),
            data=list(gene_ids))
    else:
        logging.info("Not adding \"gene_id\" to zarr output: # gene ids is 0")

    # read .npz file expression counts and add it to the expression_counts dataset
    exp_counts = np.load(args.count_matrix)
    # now convert it back to a csr_matrix object
    csr_exp_counts = sparse.csr_matrix((exp_counts['data'],
                                        exp_counts['indices'],
                                        exp_counts['indptr']),
                                       shape=exp_counts['shape'])

    if args.verbose:
        logging.info('shape of count matrix', exp_counts['shape'], csr_exp_counts.shape)
    # now create a dataset of zeros with the same dimensions as the expression count matrix
    exp_counts_group = data_group.zeros('expression',
                                        compressor=COMPRESSOR,
                                        shape=csr_exp_counts.shape,
                                        chunks=(CHUNK_ROW_SIZE, CHUNK_COL_SIZE),
                                        dtype=np.float32)

    # load a chunks from the expression count matrix and update the corresponding chunk in expression matrix in zar
    for i in range(0, csr_exp_counts.shape[0], CHUNK_ROW_SIZE):
        for j in range(0, csr_exp_counts.shape[1], CHUNK_COL_SIZE):
            # check if it is possible to make a full row chunk of data, otherwise adjust to the correct size
            p = CHUNK_ROW_SIZE
            if i + CHUNK_ROW_SIZE > csr_exp_counts.shape[0]:
                p = csr_exp_counts.shape[0] - i

            # check if it is possible to make a full row chunk of data, otherwise adjust to the correct size
            q = CHUNK_COL_SIZE
            if j + CHUNK_COL_SIZE > csr_exp_counts.shape[1]:
                q = csr_exp_counts.shape[1] - j

            # insert the chunk
            exp_counts_group[i:i + p, j:j + q] = csr_exp_counts[i:i + p, j:j + q].toarray()

    return cell_ids, gene_ids


def create_zarr_files(args):
    """This function creates the zarr file or folder structure in output_zarr_path in format file_format,
        with sample_id from the input folder analysis_output_path

    Args:
        args (argparse.Namespace): input arguments for the run
    """
    # initiate the zarr file
    root_group = init_zarr(args.sample_id, args.output_zarr_path, args.zarr_format)

    # add the expression count matrix data
    cell_ids, gene_ids = add_expression_counts(root_group['expression_matrix'], args)

    # add the the gene metrics
    add_gene_metrics(root_group['expression_matrix'], args.gene_metrics, gene_ids, args.verbose)

    # add the the cell metrics
    add_cell_metrics(root_group['expression_matrix'], args.cell_metrics, cell_ids, args.verbose)


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

    parser.add_argument('--output_path_for_zarr',
                        dest="output_zarr_path",
                        required=True,
                        help='path to .zarr file is to be created')

    parser.add_argument('--sample_id',
                        dest="sample_id",
                        default="Unknown sample",
                        help='the sample name in the bundle')

    parser.add_argument('--format',
                        dest="zarr_format",
                        default="DirectoryStore",
                        choices=["DirectoryStore", "ZipStore"],
                        help='format of the zarr file choices: [DirectoryStore, ZipStore] default: DirectoryStore')

    parser.add_argument('--verbose',
                        dest="verbose",
                        action="store_true",
                        help='whether to output verbose debugging messages')

    args = parser.parse_args()

    create_zarr_files(args)


if __name__ == '__main__':
    main()
