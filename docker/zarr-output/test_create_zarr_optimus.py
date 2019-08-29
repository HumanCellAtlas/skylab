import tempfile
import zarr
import pandas as pd

import create_zarr_optimus as target_test_code

def test_add_cell_metrics():
    path = tempfile.mkdtemp()
    store = zarr.DirectoryStore(path)
    root = zarr.group(store, overwrite=True)
    test_data_location="test/data/"
    test_metrics_file=test_data_location + "merged-cell-metrics.csv.gz"
    test_emptydrops_file=test_data_location + "empty_drops_result.csv"
    metrics_df = pd.read_csv(test_metrics_file, dtype=str)
    metrics_df = metrics_df.rename(columns={"Unnamed: 0": "cell_id"})
    sample_cell_ids=metrics_df['cell_id'][0:10]
    target_test_code.add_cell_metrics(data_group= root, metrics_file=test_metrics_file,
                                      cell_ids=sample_cell_ids, emptydrops_file=test_emptydrops_file,
                                      verbose=True);
    # Read the results back
    store_read = zarr.open(zarr.DirectoryStore(path))
    cell_metadata_uint = store_read['cell_metadata_uint'][:]
    cell_metadata_uint_name = store_read['cell_metadata_uint_name'][:]
    cell_metadata_float = store_read['cell_metadata_float'][:]
    cell_metadata_float_name = store_read['cell_metadata_float_name'][:]
    cell_metadata_bool = store_read['cell_metadata_bool'][:]
    cell_metadata_bool_name = store_read['cell_metadata_bool_name'][:]

    assert cell_metadata_uint.shape == (10,22)
    assert cell_metadata_uint_name.shape == (22,)
    assert cell_metadata_float.shape == (10,13)
    assert cell_metadata_float_name.shape == (13,)
    assert cell_metadata_bool.shape == (10,2)
    assert cell_metadata_bool_name.shape == (2,)

def run_all_tests():
    test_add_cell_metrics();

if __name__ == '__main__':
    run_all_tests();
