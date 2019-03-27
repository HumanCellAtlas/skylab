# skylab
Secondary analysis pipelines for the Human Cell Atlas.

[![GitHub Release](https://img.shields.io/github/release-pre/HumanCellAtlas/skylab.svg?label=Latest%20Release&style=flat-square&colorB=green)](https://github.com/HumanCellAtlas/skylab/releases)
[![License](https://img.shields.io/github/license/HumanCellAtlas/skylab.svg?style=flat-square)](https://github.com/HumanCellAtlas/skylab/blob/master/LICENSE)
[![Snyk Vulnerabilities for GitHub Repo (Specific Manifest)](https://img.shields.io/snyk/vulnerabilities/github/HumanCellAtlas/skylab/docker/cellranger/requirements.txt.svg?style=flat-square&label=Snyk%20Scripts%20Vulnerabilities&logo=Snyk)](https://snyk.io/test/github/HumanCellAtlas/skylab?targetFile=docker/cellranger/requirements.txt)

## Pipelines
- [cellranger](https://github.com/HumanCellAtlas/skylab/tree/master/pipelines/cellranger) secondary analysis pipeline
- [Optimus](https://github.com/HumanCellAtlas/skylab/tree/master/pipelines/optimus) secondary analysis pipeline
- [Smart-seq2](https://github.com/HumanCellAtlas/skylab/tree/master/pipelines/smartseq2_single_sample) secondary analysis pipeline

## How to run pipelines from skylab
For now, use `git clone git@github.com:HumanCellAtlas/skylab.git` and run the pipeline in Cromwell.

1. WDL and Cromwell [Documentation](https://software.broadinstitute.org/wdl/)
2. Running WDLs in [Cromwell](https://software.broadinstitute.org/wdl/documentation/execution.php)

### preemptible

Tasks on Cromwell may be run on what are known as "preemptible" machines to reduce costs by a significant amount. The catch with preemptible machines is that they may be "preempted" at any given moment--as in, google may shut down the task to re-use the resources.

Many tasks are set to automatically be `preemptible = 3`, aka they will be run on preemptible instances for up to 3 instances of preemption, after which it will be run on a non-preemptible machine. This option may be set to 0 by passing a task-level input to the workflow (i.e. Optimus.StarAlign.preemptible), causing the task to be run without using preemptible machines. This option is especially useful for long-running tasks, which have a chance to take a very long time to run due to being preempted multiple times.

### maxRetries

Some tasks have a maxRetries runtime attribute specified, with the default value set to zero. You probably shouldn't override the default and even if you do, you should do so with caution.

Setting it to an integer n greater than zero will make Cromwell retry the task up to n times if it fails for any reason. This can be useful when running a task over and over in a production setting where a high proportion of failures are due to transient problems (e.g. VM dies while job is running) that do not persist when the task is rerun. Even in that situation, it is probably best to set maxRetries to no more than 1 or 2, since if you're running in the cloud you will incur additional charges for each retry.

See the [Cromwell documentation](http://cromwell.readthedocs.io/en/develop/RuntimeAttributes/#maxretries) for more information.
