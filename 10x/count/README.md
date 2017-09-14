## 10x WDLs

These are WDLs that aim to run equivalents to `cellranger` commands.

Presently, the goal is to match results from running `cellranger` directly, but
the pipelines will likely evolve away from that over time.

### make_fastq.wdl

This WDL implements the `cellranger mkfastq` command. Its example input JSON is
example_make_fastq_input.json.

Generally, the inputs to the pipeline are an Illumina run path and a sample sheet.
The outputs are fastq files and quality metrics. The pipeline is responsible for
adjusting the sample sheet so 10x barcodes are handled correctly, running `bcl2fastq`,
collecting metrics, and merging fastqs appropriately.

Note that there are a number of parameters that may change based on different 10x
kits or experiment designs. 
