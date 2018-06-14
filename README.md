# skylab
Secondary analysis pipelines for the Human Cell Atlas.

## Pipelines
- [Smart-seq2](https://github.com/HumanCellAtlas/skylab/tree/master/pipelines/smartseq2_single_sample) secondary analysis pipeline

## How to run pipelines from skylab
For now, use `git clone git@github.com:HumanCellAtlas/skylab.git` and run the pipeline in Cromwell.

1. WDL and Cromwell [Documentation](https://software.broadinstitute.org/wdl/)
2. Running WDLs in [Cromwell](https://software.broadinstitute.org/wdl/documentation/execution.php)

### maxRetries

Some tasks have a maxRetries runtime attribute specified, with the default value set to zero. You probably shouldn't override the default and even if you do, you should do so with caution.

Setting it to an integer n greater than zero will make Cromwell retry the task up to n times if it fails for any reason. This can be useful when running a task over and over in a production setting where a high proportion of failures are due to transient problems (e.g. VM dies while job is running) that do not persist when the task is rerun. Even in that situation, it is probably best to set maxRetries to no more than 1 or 2, since if you're running in the cloud you will incur additional charges for each retry.

See the [Cromwell documentation](http://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#maxretries) for more information.
