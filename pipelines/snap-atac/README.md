# SnapATAC

## Requirements

- `cromwell-tools`
- Service account key
- For WDL validation, the bash scripts expect `womtool-40.jar` to be in `~/Applications/` (fixme)

## How to Run Examples

### Human

Below will submit a job to CaaS:

```bash
$ ./snap-atac.submit.sh \
    -i snap-atac.inputs.human.json \
    -k ~/pipelines_surge_credentials.json
```

### Mouse

Below will submit a job to CaaS:

```bash
$ ./snap-atac.submit.sh \
    -i snap-atac.inputs.mouse.json \
    -k ~/pipelines_surge_credentials.json
```

### Barcode Attachment

SnapATAC only takes barcoded-attached FASTQ as input. Run this if you have four FASTQ (I5, I7, R1, R2) as input:

```bash
$ bc-attacher.submit.sh \
    -k ~/pipelines_surge_credentials.json
```

which will generate two R1, R2 FASTQ fiels.

## WDL Validation

Below will validate your WDLs:

```bash
$ ./snap-atac.validate.sh
```

```bash
$ ./bc-attacher.validate.sh
```
