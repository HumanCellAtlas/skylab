# skylab
Secondary analysis pipelines for the Human Cell Atlas.

## Pipelines
- [Smart-seq2](https://github.com/HumanCellAtlas/skylab/tree/master/smartseq2_single_sample) secondary analysis pipeline

**Bold indicates that this pipeline is blessed**

## How to run pipelines from skylab
For now, use `git clone git@github.com:HumanCellAtlas/skylab.git` and run the pipeline in Cromwell.

1. WDL and Cromwell [Documentation](https://software.broadinstitute.org/wdl/)
2. Running WDLs in [Cromwell](https://software.broadinstitute.org/wdl/documentation/execution.php)

**[TODO] Update this section with better instructions on how to easily try our pipelines.**

## [TODO] Blessed HCA Pipelines in Methods Repositories
- Agora
- Dockstore

## Example on how to add a pipeline (workflow) to Dockstore and run it
After Dockstore registration has been completed, register the pipeline (workflow) on Dockstore with the following fields:
  - "HumanCellAtlas/skylab" as the "Source Code Repository"
  - "/pipelines/smartseq2_single_sample/SmartSeq2SingleSample.wdl" as the "Workflow Path"
  - "/pipelines/smartseq2_single_sample/dockstore_SmartSeq2SingleSampleExample.json" as the "Test Parameter File Path"
After registering, publish the pipeline (workflow) and download the dockstore_SmartSeq2SingleSampleExample.json if it's not present already. 
It can then be ran using the [Dockstore CLI](https://dockstore.org/quick-start):

`dockstore workflow launch --entry github.com/HumanCellAtlas/skylab:master --json dockstore_SmartSeq2SingleSampleExample.json`

## [TODO] Tagged releases
Pipelines are released in tagged releases to track changes in version of pipeline, and to track when a pipeline is blessed.

## [TODO] Automated testing options
## [TODO] Onboarding pipelines
