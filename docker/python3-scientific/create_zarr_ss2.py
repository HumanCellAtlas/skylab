import argparse
import csv
import numpy
import zarr
from zarr import Blosc


COMPRESSOR = Blosc(cname='lz4', clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)

def init_zarr(sample_id, path, file_format):
    """Initializes the zarr output
    Args:
        sample_id (str): sample or cell id
        path (str): path to the zarr output
        fileformat (str): zarr file format [ DirectoryStore, ZipStore]
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
    root.attrs['README'] = ("The schema adopted in this zarr store may undergo "
                            "changes in the future")
    root.attrs['sample_id'] = sample_id

    return root


def read_and_convert_qc(data_group, qc_paths):
    """Converts the QC of Smart Seq2 gene file pipeline outputs to zarr file
    Args:
        data_group (zarr.hierarchy.Group): zarr group to add the datasets
        qc_path (str): path to the QCs csv
    """
    # read the QC values
    qc_path = [p for p in qc_paths if p.endswith("_QCs.csv")][0]

    with open(qc_path, 'r') as f:
        qc_values = [row for row in csv.reader(f)]

    metadata_labels = qc_values[0][1:]
    metadata_values = qc_values[2][1:]
    cell_id = qc_values[2][0]

    string_metadata = {}
    numeric_metadata = {}

    for label, value in zip(metadata_labels, metadata_values):

        # See if this is numeric or string
        numeric_value = None
        try:
            numeric_value = float(value)
        except ValueError:
            try:
                numeric_value = float(value.strip("%"))/100
            except ValueError:
                pass

        if numeric_value is not None:
            numeric_metadata[label] = numeric_value
        else:
            string_metadata[label] = value

    # Metrics
    # Write the string and numeric metadata separately
    sorted_string_labels = sorted(string_metadata.keys())
    sorted_string_values = [string_metadata[m] for m in sorted_string_labels]
    data_group.create_dataset(
        "cell_metadata_string_name",
        shape=(len(string_metadata),),
        compressor=COMPRESSOR,
        dtype='<U40',
        chunks=(len(string_metadata),),
        data=sorted_string_labels
    )
    data_group.create_dataset(
        "cell_metadata_string",
        shape=(1, len(string_metadata)),
        compressor=COMPRESSOR,
        dtype='<U40',
        chunks=(1, len(string_metadata)),
        data=[sorted_string_values]
    )

    sorted_numeric_labels = sorted(numeric_metadata.keys())
    sorted_numeric_values = [numeric_metadata[m] for m in sorted_numeric_labels]
    data_group.create_dataset(
        "cell_metadata_numeric_name",
        shape=(len(numeric_metadata),),
        compressor=COMPRESSOR,
        dtype="<U40",
        chunks=(len(numeric_metadata),),
        data=sorted_numeric_labels
    )
    data_group.create_dataset(
        "cell_metadata_numeric",
        shape=(1, len(numeric_metadata)),
        compressor=COMPRESSOR,
        dtype="f4",
        chunks=(1, len(numeric_metadata)),
        data=[sorted_numeric_values]
    )

    # Cell IDs
    data_group.create_dataset(
        "cell_id",
        shape=(1,),
        compressor=COMPRESSOR,
        dtype="<U40",
        chunks=(1,),
        data=[cell_id]
    )


def read_and_convert_expression(data_group, input_path):
    """Converts the Smart Seq2 gene file pipeline outputs to zarr file
    Args:
        data_group (zarr.hierarchy.Group): databroup object for the zarr
        input_path (str): file where the SS2 pipeline expression counts are
    """
    reader = csv.DictReader(open(input_path), delimiter="\t")

    expression_values = {}
    for row in reader:
        expression_values[row["gene_id"]] = float(row["TPM"])

    sorted_gene_ids = sorted(expression_values.keys())
    sorted_tpms = [expression_values[g] for g in sorted_gene_ids]

    # TPM
    data_group.create_dataset(
        "expression",
        shape=(1, len(expression_values)),
        compressor=COMPRESSOR,
        dtype=numpy.float32,
        chunks=(1, len(expression_values)),
        data=[sorted_tpms]
    )

    # Gene IDs
    data_group.create_dataset(
        "gene_id",
        shape=(len(expression_values),),
        compressor=COMPRESSOR,
        dtype="<U40",
        chunks=(len(expression_values),),
        data=sorted_gene_ids
    )


def create_zarr_files(sample_id, qc_files, rsem_genes_results_file,
                      output_zarr_path, file_format="DirectoryStore"):
    """This function creates the zarr file or folder structure in output_zarr_path in
       format file_format, with sample_id from the input folder analysis_output_path
    Args:
        sample_id (str): sample or cell id
        qc_analysis_output_files_string (str): a string with the file names in the QCGroup of SS2
            pipeline output, separated by commas
        rsem_genes_results_file (str): the file for the expression count
        output_zarr_path (str): location of the output zarr
        file_format (str): zarr file format  [DirectoryStore or ZipStore ] Default: DirectoryStore
    """
    # initiate the zarr file
    root_group = init_zarr(sample_id, output_zarr_path, file_format)

    # add the the gene/TMP counts
    read_and_convert_expression(root_group, rsem_genes_results_file)

    # add the the QC metrics
    read_and_convert_qc(root_group, qc_files)


def main():
    description = """This script converts the some of the SmartSeq2 pipeline outputs in to
                   zarr format (https://zarr.readthedocs.io/en/stable/) relevant output. 
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--qc_files',
                        nargs="+",
                        help=('the grouped QC files from the GroupQCOutputs task of SS2 '
                              'Single Sample workflow'))

    parser.add_argument('--rsem_genes_results',
                        help='path to the folder containing the files to be added to the zarr')

    parser.add_argument('--output_zarr_path',
                        help='path where the zarr file is to be created')

    parser.add_argument('--sample_id',
                        default="Unknown sample",
                        help='the sample name in the bundle')

    parser.add_argument('--format',
                        default="DirectoryStore",
                        help=('format of the zarr file choices: '
                              '[DirectoryStore, ZipStore] default: DirectoryStore'))
    args = parser.parse_args()

    create_zarr_files(args.sample_id, args.qc_files, args.rsem_genes_results,
                      args.output_zarr_path, args.format)


if __name__ == '__main__':
    main()
