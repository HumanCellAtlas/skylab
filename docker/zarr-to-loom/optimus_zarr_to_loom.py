#!/usr/bin/env python3
import os
import sys

import zarr
import loompy
import scipy as sc
from scipy.sparse import coo_matrix
import argparse
import numpy as np


# Custom Exception
class UnexpectedInputFormat(Exception):
    pass


def main():
    # Parse the arguments
    description = """This script converts ZARR optimus output to loom format"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-zarr', dest="input_zarr_path", required=True, help="Path to input ZARR file")
    parser.add_argument('--output-loom', dest="output_loom_path", required=True, help="Path to output loom file")
    parser.add_argument('--sample-id', dest="sample_id", required=True, help="Sample identifier")
    args = parser.parse_args()

    input_zarr_path = args.input_zarr_path
    output_loom_path = args.output_loom_path
    sample_id = args.sample_id

    # Checks on inputs
    if not os.path.isdir(input_zarr_path):
        sys.exit("Error: the input ZARR path is not a directory.")
    if os.path.exists(output_loom_path):
        sys.exit("Error: The output loom file exists!")

    # Open the ZARR
    store = zarr.DirectoryStore(input_zarr_path)
    root = zarr.open(store)

    # Get the expression matrix
    # expression matrix in numpy ndarray format (dense)
    # NOTE: If memory is limiting this could be done by chunk
    nrows = root[f"/{sample_id}.zarr/expression"].shape[0]
    ncols = root[f"/{sample_id}.zarr/expression"].shape[1]

    expr_sp = sc.sparse.coo_matrix((nrows, ncols), np.float32)

    iter_row_count = 100;

    xcoord = []
    ycoord = []
    value = []

    chunk_row_size = 10000
    chunk_col_size = 10000

    for i in range(0, nrows, chunk_row_size):
        for j in range(0, ncols, chunk_col_size):
            p = chunk_row_size
            if i + chunk_row_size > nrows:
                p = nrows - 1
            q = chunk_col_size
            if j + chunk_col_size > ncols:
                q = ncols - j
            expr_sub_row_coo = sc.sparse.coo_matrix(root[f"/{sample_id}.zarr/expression"][i:(i+p), j:(j+q)])
            for k in range(0, expr_sub_row_coo.data.shape[0]):
                xcoord.append(expr_sub_row_coo.row[k] + i)
                ycoord.append(expr_sub_row_coo.col[k] + j)
                value.append(expr_sub_row_coo.data[k])

    xcoord = np.asarray(xcoord)
    ycoord = np.asarray(ycoord)
    value = np.asarray(value)

    expr_sp_t = sc.sparse.coo_matrix((value, (ycoord, xcoord)), shape=(expr_sp.shape[1], expr_sp.shape[0]))

    del xcoord
    del ycoord
    del value

    # ROW/GENE Metadata

    # Check that the first gene metadata column is the gene name as expected
    if not (root[f"/{sample_id}.zarr/gene_metadata_string_name"][:] == 'gene_name')[0]:
        raise UnexpectedInputFormat("The first gene metadata item is not the gene_name");

    # Prepare row attributes (gene attributes)
    # Follow loom file Conventions
    row_attrs = {
        "Gene": root[f"/{sample_id}.zarr/gene_metadata_string"][:][0,],
        "Accession": root[f"/{sample_id}.zarr/gene_id"][:]}

    numeric_field_names = root[f"/{sample_id}.zarr/gene_metadata_numeric_name"][:]

    # Generate with a list comprehension
    def add_to_gene_meta_by_index(i):
        name = numeric_field_names[i]
        data = root[f"/{sample_id}.zarr/gene_metadata_numeric"][:][:, i]
        row_attrs[name] = data

    [add_to_gene_meta_by_index(i) for i in range(0, numeric_field_names.shape[0])]

    # COLUMN/CELL Metadata
    col_attrs = dict()
    col_attrs["CellID"] = root[f"/{sample_id}.zarr/cell_id"][:]
    bool_field_names = root[f"/{sample_id}.zarr/cell_metadata_bool_name"][:]

    def add_to_cell_meta_bool_by_index(i):
        name = bool_field_names[i]
        data = root[f"/{sample_id}.zarr/cell_metadata_bool"][:][:, i]
        col_attrs[name] = data

    [add_to_cell_meta_bool_by_index(i) for i in range(0, bool_field_names.shape[0])]

    float_field_names = root[f"/{sample_id}.zarr/cell_metadata_float_name"][:]

    def add_to_cell_meta_float_by_index(i):
        name = float_field_names[i]
        data = root[f"/{sample_id}.zarr/cell_metadata_float"][:][:, i]
        col_attrs[name] = data

    [add_to_cell_meta_float_by_index(i) for i in range(0, float_field_names.shape[0])]

    # Generate the loom file
    loompy.create(output_loom_path, expr_sp_t, row_attrs, col_attrs)


if __name__ == '__main__':
    main()
