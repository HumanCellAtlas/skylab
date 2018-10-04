import argparse
import csv
import glob
import json
import os, sys
import re, logging

import crimson.picard
import numcodecs.blosc
import numpy
import zarr

""" key,  dataset name,   suffix of file,  metadata description 
    this is a list of groups to add to the zarr file, which can be 
    further extended
"""
_FILES = {
  "QCs" : [ "QCs", "*_QCs.csv",  "QC data" ], 
  "genes_results"  : [ "TPM", "*_rsem.genes.results", "count for the genes and also isoforms" ],
}

COMPRESSOR = numcodecs.blosc.Blosc(cname='lz4', clevel=5, shuffle=1, blocksize=0)

def init_zarr(sample_id, path, fileformat):
    """  Creates an empty zarr file  in the path
        Args:
            path:  path to the zarr 
        Returns:
            
        Raises:
    """
    if fileformat=="DirectoryStore":
        store = zarr.DirectoryStore(path)

    if fileformat=="ZipStore":
       store= zarr.ZipStore(path, mode='w')

    # create the root group
    root = zarr.group(store, overwrite=True)

    # add some readme for the user
    root.attrs['README'] = "The schema adopted in this zarr store may undergo changes in future"
    root.attrs['sample_id'] = sample_id

    # now iterate through list of expected groups and create them
    for data in _FILES:
      datagroup = root.create_group(data, overwrite=True)
      datagroup.attrs['description'] = _FILES[data][2]

    return root


def read_and_convert_QC(root, input_files_string):
    """Converts the QC of Smart Seq2 gene file pipeline outputs to zarr file
        Args:
            root: root object for the zarr
            analysis_output_path: path where the SS2 pipeline are found 
        Returns:
            True if success else False
        Raises:

    """
    
    # split the files string by , 
    # crate an array of the file names containing the suffix _QCs.csv 
    relevant_qc_files =  [  x.strip() for x in input_files_string.split(',')  if re.search(r'_QCs.csv$', x) ]

    # we expect only one, write an warning 
    if len(relevant_qc_files) !=1:
       logging.warning(sys.stderr.write("WARNING: Multiple files ending with _QCs.csv. Expected only one"))


    # read the QC values
    QC_values = []
    with open(relevant_qc_files[0], 'r') as f:
       reader = csv.reader(f)
       for row in reader:
          QC_values.append(row)

    datagroup = root['QCs']

    # Metrics
    """ The first value is the row name so we remove it"""
    
    metrics_array=datagroup.create_dataset(
        "QCs",
        shape=(len(QC_values[0][1:]), ),
        dtype='<U40',
        chunks=(len(QC_values[0][1:]), ),
        data = QC_values[0][1:]
     )

    # Class
    metrics_array=datagroup.create_dataset(
        "Class",
        shape=(len(QC_values[0][1:]), ),
        dtype='<U40',
        chunks=(len(QC_values[0][1:]), ),
        data = QC_values[0][1:] 
     )

   # Values, which are the values of the metrics
    metrics_array=datagroup.create_dataset(
        "Values",
        shape=(len(QC_values[2][1:]), ),
        dtype='<U40',
        chunks=(len(QC_values[2][1:]), ),
        data = QC_values[2][1:]
     )
      

def read_and_convert_expression(root, input_path):
    """Converts the Smart Seq2 gene file pipeline outputs to zarr file
        Args:
            root: root object for the zarr
            analysis_output_path: path where the SS2 pipeline are found 
        Returns:
            True if success else False
        Raises:
            FileNotFoundError
    """
    #expression_path = glob.glob(os.path.join(input_path, _FILES['genes_results'][1]))[0]
    reader = csv.DictReader(open(input_path), delimiter="\t")

    expression_values = {}
    for row in reader:
        expression_values[row["gene_id"]] = float(row["TPM"])

    datagroup = root['genes_results']

    # Gene IDs
    gene_id_array = datagroup.create_dataset(
        "gene_id",
        shape=(len(expression_values),),
        dtype="<U40",
        chunks=(len(expression_values), ), 
        data=list(expression_values.keys())
    )

    # TPM
    expression_array=datagroup.create_dataset(
        "TPM",
        shape=(len(expression_values), ),
        dtype=numpy.float32,
        chunks=(len(expression_values), ),
        data= list(expression_values.values())
     )
      


def create_zarr_files(sample_id,qc_analysis_output_files_string, rsem_genes_results_file,  output_zarr_path, fileformat):
    """ This function creates the zarr file in output_zarr_path in format fileformat,
        with sample_id from the input folder analysis_output_path  """

    #initiate the zarr file 
    root = init_zarr(sample_id, output_zarr_path, fileformat)

    # add the the gene/TMP counts
    read_and_convert_expression(root, rsem_genes_results_file)

    # add the the QC metrics
    read_and_convert_QC(root, qc_analysis_output_files_string)

def main():
    description = """This script converts the some of the SmartSeq2 pipeline outputs in to
                   zarr format (https://zarr.readthedocs.io/en/stable/) relevant output. 
                   This script can be used as a module or run as a command line script."""

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--qc_analysis_output_files_string', dest="qc_analysis_output_files_string", 
                        help='a string separated by commas with a list of files from the GroupQCOutputs step of SS2 pipeline'
                       )

    parser.add_argument('--rsem_genes_results', dest="rsem_genes_results", 
                        help='path to the folder containing the files to be added to the zarr'
                       )
    parser.add_argument('--output_path_for_zarr', dest="output_zarr_path",
                        help='path where the zarr file is to be created'
                       )
    parser.add_argument('--sample_id', dest="sample_id", default="Unknown sample",
                        help='the sample name in the bundle'
                       )

    parser.add_argument('--format', dest="format", default="DirectoryStore",
                        help='format of the zarr file choices: [DirectoryStore, ZipStore] default: DirectoryStore'
                       )
    args = parser.parse_args()

    create_zarr_files(args.sample_id, args.qc_analysis_output_files_string, args.rsem_genes_results, args.output_zarr_path, args.format)
    #add_cell_metadata(root, args.analysis_output_path)
    #add_gene_metadata(root, args.analysis_output_path)

if __name__ == '__main__':
    main()
