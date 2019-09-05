import argparse
import csv
import gzip
import re
import numpy as np
import zarr
from scipy import sparse
from zarr import Blosc
import pandas as pd
import logging


COMPRESSOR = Blosc(cname='lz4', clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)

# the number of rows in a chunk for expression counts
CHUNK_ROW_SIZE = 10000
CHUNK_COL_SIZE = 10000
logging.basicConfig(level=logging.INFO)


def init_zarr(sample_id, path, file_format, schema_version):
    """Initializes the zarr output.

    Args:
        sample_id (str): sample or cell id
        path (str): path to the zarr output
        file_format (str): zarr file format [DirectoryStore, ZipStore]
        schema_version (str): version string of this output to allow for parsing of future changes

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

    #root.attrs['README'] = "The schema adopted in this zarr store may undergo changes in the future"
    root.attrs['sample_id'] = sample_id
    root.attrs['optimus_output_schema_version'] = schema_version

    # Create the expression_matrix group
    #root.create_group("expression_matrix", overwrite=True);

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


def add_cell_metrics(data_group, metrics_file, cell_ids, emptydrops_file, verbose=False,):
    """Converts cell metrics from the Optimus pipeline to zarr file

    Args:
        data_group (zarr.hierarchy.Group): datagroup object for the zarr
        input_path (str): file containing gene metrics name and values
        cell_ids (list): list of cell ids
        verbose (bool): whether to output verbose messages for debugging purposes
        emptydrops_path (str): emptydrops csv file
    """

    # Read the csv input files
    metrics_df = pd.read_csv(metrics_file, dtype=str)
    emptydrops_df = pd.read_csv(emptydrops_file, dtype=str)

    # Check that input is valid
    if (metrics_df.shape[0] == 0 or metrics_df.shape[1] == 0):
        logging.error("Cell metrics table is not valid")
        raise ValueError()
    if (emptydrops_df.shape[0] == 0 or emptydrops_df.shape[1] == 0):
        logging.error("EmptyDrops table is not valid")
        raise ValueError()

    # Rename cell columns for both datasets to cell_id
    emptydrops_df = emptydrops_df.rename(columns={"CellId": "cell_id"})
    metrics_df = metrics_df.rename(columns={"Unnamed: 0": "cell_id"})

    # Drop first row that contains non-cell information from metrics file, this contains aggregate information
    metrics_df = metrics_df.iloc[1:]

    # Prefix emptydrops column names (except the key cell_id)
    colnames = list(emptydrops_df.columns)
    newcolnames = ["emptydrops_" + s for s in colnames]
    namemap = dict(zip(colnames, newcolnames))
    # Do not map the cell_id as it will be used for the merge
    del namemap["cell_id"]
    emptydrops_df = emptydrops_df.rename(columns=namemap)

    # Confirm that the emptydrops table is a subset of the cell metadata table, fail if not
    if (not emptydrops_df.cell_id.isin(metrics_df.cell_id).all()):
        logging.error("Not all emptydrops cells can be found in the metrics table.")
        raise Exception("Not all emptydrops cells can be found in the metrics table.")

    # Merge the two tables
    merged_df = metrics_df.merge(emptydrops_df, on="cell_id", how="outer")

    # Order the cells by merging with cell_ids
    cellorder_df = pd.DataFrame(data={'cell_id': cell_ids})
    final_df = cellorder_df.merge(merged_df, on='cell_id', how="left")

    # Split the pandas DataFrame into different data types for storing in the ZARR
    FloatColumnNames = [# UInt
                        "n_reads", "noise_reads", "perfect_molecule_barcodes",
                        "reads_mapped_exonic", "reads_mapped_intronic", "reads_mapped_utr",
                        "reads_mapped_uniquely", "reads_mapped_multiple", "duplicate_reads",
                        "spliced_reads", "antisense_reads", "n_molecules", "n_fragments",
                        "fragments_with_single_read_evidence", "molecules_with_single_read_evidence",
                        "perfect_cell_barcodes", "reads_mapped_intergenic",
                        "reads_unmapped", "reads_mapped_too_many_loci",
                        "n_genes", "genes_detected_multiple_observations",
                        "emptydrops_Total",
                        # Float32
                        "molecule_barcode_fraction_bases_above_30_mean",
                        "molecule_barcode_fraction_bases_above_30_variance",
                        "genomic_reads_fraction_bases_quality_above_30_mean",
                        "genomic_reads_fraction_bases_quality_above_30_variance",
                        "genomic_read_quality_mean",
                        "genomic_read_quality_variance",
                        "reads_per_fragment",
                        "fragments_per_molecule",
                        "cell_barcode_fraction_bases_above_30_mean",
                        "cell_barcode_fraction_bases_above_30_variance",
                        "emptydrops_LogProb",
                        "emptydrops_PValue",
                        "emptydrops_FDR"]
    BoolColumnNames = ["emptydrops_Limited", "emptydrops_IsCell"]

    # Split the dataframe
    final_df_float = final_df[FloatColumnNames]
    final_df_bool = final_df[BoolColumnNames]

    # Data types for storage
    header_datatype = "<U40"  # little-endian 40 char unicode
    float_store_datatype = np.float32  # machine independent 32 bit float
    bool_store_datatype = np.bool  # boolean

    # Do format conversions
    final_df_float = final_df_float.apply(pd.to_numeric)

    # Create metadata tables and their headers for float
    data_group.create_dataset("cell_metadata_float_name", shape=[final_df_float.shape[1], 1],
                              compressor=COMPRESSOR, dtype=header_datatype, data=final_df_float.columns.astype(str))
    data_group.create_dataset("cell_metadata_float", shape=final_df_float.shape, compressor=COMPRESSOR,
                              dtype=float_store_datatype, data=final_df_float.to_numpy(dtype=float_store_datatype))
    if verbose:
        logging.info("Added cell metadata_float with {} rows and {} columns".format(
            final_df_float.shape[0], final_df_float.shape[1]))

    # Create metadata tables and their headers for bool
    data_group.create_dataset("cell_metadata_bool_name", shape=[final_df_bool.shape[1], 1],
                              compressor=COMPRESSOR, dtype=header_datatype, data=final_df_bool.columns.astype(str))
    data_group.create_dataset("cell_metadata_bool", shape=final_df_bool.shape, compressor=COMPRESSOR,
                              dtype=bool_store_datatype, data=final_df_bool.to_numpy(dtype=bool_store_datatype))
    if verbose:
        logging.info("Added cell metadata_bool with {} rows and {} columns".format(
            final_df_bool.shape[0], final_df_bool.shape[1]))

def create_gene_id_name_map(gtf_file):
    """ Creates a map from gene_id to gene_name by reading in the GTF file

    Args:
        gtf_file (str): annotation file

    Return:
        gene_id_name_map (Dict[str, str]): dictonary gene ids to gene names
    """
    gene_id_name_map = {}

    # loop through the lines and find the gene_id and gene_name pairs
    with gzip.open(gtf_file, 'rt') if gtf_file.endswith(".gz") else open(gtf_file, 'r') as fpin:
        for _line in fpin:
            line = _line.strip()
            gene_id_res = re.search(r'gene_id ([^;]*);', line)
            gene_name_res = re.search(r'gene_name ([^;]*);', line)
    
            if gene_id_res and gene_name_res:
                gene_id = gene_id_res.group(1).replace('"', '')
                gene_name = gene_name_res.group(1).replace('"', '')
                gene_id_name_map[gene_id] = gene_name

    return gene_id_name_map


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

    # add the gene names for the gene ids 
    if args.annotation_file and len(gene_ids):
        gene_id_name_map = create_gene_id_name_map(args.annotation_file)
        gene_names = [gene_id_name_map.get(gene_id, "") for gene_id in gene_ids]

        # insert the name of the "gene id to gene name map"  gene string metadata
        data_group.create_dataset('gene_metadata_string_name',
                                  compression=COMPRESSOR,
                                  dtype='<U40',
                                  chunks=(1, ),
                                  data=["gene_name"])

        # insert the array of  gene name as an implicit "gene id to gene name map"  to the gene metadata
        data_group.create_dataset('gene_metadata_string',
                                  shape=(1, len(gene_ids),),
                                  compression=COMPRESSOR,
                                  dtype='<U40',
                                  chunks=(1, len(gene_ids)),
                                  data=[gene_names])

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
    version = "1.0.0"

    # initiate the zarr file
    root_group = init_zarr(args.sample_id, args.output_zarr_path, args.zarr_format, schema_version=version)

    # add the expression count matrix data
    cell_ids, gene_ids = add_expression_counts(root_group, args)

    # add the the gene metrics
    add_gene_metrics(root_group, args.gene_metrics, gene_ids, args.verbose)

    # add the the cell metrics
    add_cell_metrics(root_group, args.cell_metrics, cell_ids,
                     args.empty_drops_file, args.verbose)

def main():
    description = """This script converts the some of the Optimus outputs in to
                   zarr format (https://zarr.readthedocs.io/en/stable/) relevant output. 
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)

    parser.add_argument('--empty_drops_file',
                        dest="empty_drops_file",
                        required=True,
                        help="A csv file with the output of the emptyDrops step in Optimus")

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

    parser.add_argument('--annotation_file',
                        dest="annotation_file",
                        default=None,
                        required=False,
                        help='annotation file in GTF format')

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
