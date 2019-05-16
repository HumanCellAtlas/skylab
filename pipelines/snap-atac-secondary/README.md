# SnapATAC Secondary Analysis

This takes a `.snap` file, runs the secondary analysis, and produces a tSNE, UMAP, and etc in PDF.

## Requirements

- `cromwell-tools`
- Service account key
- For WDL validation, the bash scripts expect `womtool-40.jar` to be in `~/Applications/` (fixme)

## How to Run Examples

The example is 10x fresh cortex from adult mouse brain (P50):
https://support.10xgenomics.com/single-cell-atac/datasets/1.1.0/atac_v1_adult_brain_fresh_5k

Below will submit a job to CaaS:

```bash
$ ./submit.sh \
    -k ~/pipelines_surge_credentials.json
```

## WDL Validation

```bash
$ ./validate.sh
```
