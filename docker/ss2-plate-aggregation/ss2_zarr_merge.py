#!/usr/bin/env python3

import os
import zarr
import argparse
import numpy as np
from zarr import Blosc

COMPRESSOR = Blosc(cname='lz4', clevel=5, shuffle=Blosc.SHUFFLE, blocksize=0)


class MismatchingInputHeader(Exception):
    pass


class PathNotDirectory(Exception):
    pass


def main():
    description = """Merge the outputs of multiple SS2 pipeline runs into a single Zarr file"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--input-zarr-directory',
                        dest='input_zarr_dir',
                        required=True,
                        help="Path to input zarr directory")
    parser.add_argument('--output-zarr-file',
                        dest='output_zarr_file',
                        required=True,
                        help="Path to output zarr file")
    parser.add_argument('--plate-sample-id',
                        dest='plate_sample_id',
                        required=True,
                        help="Plate sample id")
    parser.add_argument('--check-all-headers',
                        dest='check_all_headers',
                        action='store_true',
                        help="Check that the input from all files match")
    args = parser.parse_args()

    # The list of ZARR files that we need to merge
    zarr_file_list = os.listdir(args.input_zarr_dir)

    first_cell = True
    number_of_cells = len(zarr_file_list)
    current_cell_count = 0

    for zarr_cell_dir in zarr_file_list:
        abs_path_dir = os.path.join(args.input_zarr_dir, zarr_cell_dir)
        if not os.path.isdir(abs_path_dir):
            raise PathNotDirectory("Path in zarr directory is not a directory")
        store = zarr.DirectoryStore(abs_path_dir)
        root = zarr.open(store)
        if first_cell:
            output_store = zarr.DirectoryStore(args.output_zarr_file)
            output_root = zarr.group(output_store, overwrite=True)
            output_root.attrs['README'] = ("The schema adopted in this zarr store may undergo "
                                           "changes in the future")
            output_root.attrs['sample_id'] = args.plate_sample_id
            # Get the number of genes we have to store
            number_of_genes = len(root.gene_id[:])
            number_of_numeric_metadata = len(root.cell_metadata_numeric_name[:])
            number_of_string_metadatata = len(root.cell_metadata_string_name[:])
            # Create the cell_id dataset
            group_cell_id = output_root.create_dataset(
                "cell_id",
                shape=(number_of_cells,),
                compressor=COMPRESSOR,
                dtype="<U40",
                chunks=(number_of_cells,)
            )
            # Create the gene_id dataset
            group_gene_id = output_root.create_dataset(
                "gene_id",
                shape=(number_of_genes,),
                compressor=COMPRESSOR,
                dtype="<U40",
                chunks=(number_of_genes,),
                data=root.gene_id[:]
            )
            # Create the expression value store
            group_expression = output_root.create_dataset(
                "expression",
                shape=(number_of_cells, number_of_genes),
                compressor=COMPRESSOR,
                dtype=np.float32,
                chunks=(1, number_of_genes)
            )
            group_numeric_metadata_name = output_root.create_dataset(
                "cell_metadata_numeric_name",
                shape=(number_of_numeric_metadata,),
                compressor=COMPRESSOR,
                dtype="<U40",
                chunks=(number_of_numeric_metadata,),
                data=root.cell_metadata_numeric_name[:]
            )
            group_string_metadata_name = output_root.create_dataset(
                "cell_metadata_string_name",
                shape=(number_of_string_metadatata,),
                compressor=COMPRESSOR,
                dtype="<U40",
                chunks=(number_of_string_metadatata,),
                data=root.cell_metadata_string_name[:]
            )
            group_numeric_metadata = output_root.create_dataset(
                "cell_metadata_numeric",
                shape=(number_of_cells, number_of_numeric_metadata),
                compressor=COMPRESSOR,
                dtype="f4",
                chunks=(1, number_of_numeric_metadata)
            )
            group_string_metadata = output_root.create_dataset(
                "cell_metadata_string",
                shape=(number_of_cells, number_of_string_metadatata),
                compressor=COMPRESSOR,
                dtype="<U40",
                chunks=(1, number_of_string_metadatata)
            )
            first_cell = False
        else:
            if (args.check_all_headers):
                    if not np.array_equal(group_gene_id[:], root.gene_id[:]):
                        raise MismatchingInputHeader("Gene ids didn't match")
                    if not np.array_equal(group_numeric_metadata_name[:], root.cell_metadata_numeric_name[:]):
                        raise MismatchingInputHeader("Numeric metadata names didn't match")
                    if not np.array_equal(group_string_metadata_name[:], root.cell_metadata_string_name[:]):
                        raise MismatchingInputHeader("String metadata names didn't match")
        # Save the cell name
        group_cell_id[current_cell_count,] = root.cell_id[0]
        group_expression[current_cell_count, :] = root.expression[0, :]
        group_numeric_metadata[current_cell_count, :] = root.cell_metadata_numeric[0, :]
        group_string_metadata[current_cell_count, :] = root.cell_metadata_string[0, :]
        # Increment cell count
        current_cell_count += 1


if __name__ == '__main__':
    main()