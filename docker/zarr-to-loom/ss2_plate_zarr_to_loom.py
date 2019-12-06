#!/usr/bin/env  python3


import os
import sys
import zarr
import loompy
import scipy as sc
from scipy.sparse import coo_matrix
import argparse
import numpy as np

def main():
    description = """This script converts SS2 pipeline zarr output into loom"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-zarr', dest='input_zarr_path', required=True, help="Path to input zarr file")
    parser.add_argument('--output-loom', dest='output_loom_path', required=True, help="Path to output loom file")
    parser.add_argument('--sample-id', dest='sample_id', required=True, help="Sample identifier")
    args = parser.parse_args()

    if not os.path.isdir(args.input_zarr_path):
        sys.exit("Error: the input zarr path is not a directoyr")
    if os.path.exists(args.output_loom_path):
        sys.exit("Error: the output loom file exists")

    store = zarr.DirectoryStore(args.input_zarr_path)
    root = zarr.open(store)

    nrows = root["/output_zarr/expression"].shape[0]
    ncols = root["/output_zarr/expression"].shape[1]
    expr_sp = sc.sparse.coo_matrix((nrows, ncols), np.float32)

    xcoord = []
    ycoord = []
    value = []

    expr_coo = sc.sparse.coo_matrix(root["/output_zarr/expression"][:])
    for k in range(0, expr_coo.data.shape[0]):
        xcoord.append(expr_coo.row[k])
        ycoord.append(expr_coo.col[k])
        value.append(expr_coo.data[k])

    xcoord = np.asarray(xcoord)
    ycoord = np.asarray(ycoord)
    value = np.asarray(value)

    expr_sp_t = sc.sparse.coo_matrix((value, (ycoord, xcoord)), shape=(expr_sp.shape[1], expr_sp.shape[0]))

    del xcoord
    del ycoord
    del value

    row_attrs = {
        "Gene": root.output_zarr.gene_id[:]
    }

    col_attrs = dict()
    col_attrs["CellID"] = root.output_zarr.cell_id[:]


    loompy.create(args.output_loom_path, expr_sp_t, row_attrs, col_attrs)


if __name__ == '__main__':
    main()