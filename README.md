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

After registering, publish the pipeline (workflow) and it will be available via the [Dockstore's GA4GH TRS endpoints](https://dockstore.org:8443/static/swagger-ui/index.html) and [Dockstore CLI](https://dockstore.org/quick-start).  

### Dockstore CLI
Download the dockstore_SmartSeq2SingleSampleExample.json and run the pipeline (workflow) using:
```
$ wget https://dockstore.org:8443/api/ga4gh/v2/tools/%23workflow%2Fgithub.com%2FHumanCellAtlas%2Fskylab/versions/master/PLAIN_WDL/tests
$ dockstore workflow launch --entry github.com/HumanCellAtlas/skylab:master --json tests
```

### TRS and Cromwell
See the list of files available using GET /api/ga4gh/v2/tools/{id}/versions/{version_id}/{type}/files:
```
$ curl -X GET "https://dockstore.org:8443/api/ga4gh/v2/tools/%23workflow%2Fgithub.com%2FHumanCellAtlas%2Fskylab/versions/master/WDL/files" -H  "accept: application/json"
```
This will indicate which files are needed to run the pipeline (workflow).  It should look like:
```
[
  {
    "file_type": "PRIMARY_DESCRIPTOR",
    "path": "/pipelines/smartseq2_single_sample/SmartSeq2SingleSample.wdl"
  },
  {
    "file_type": "TEST_FILE",
    "path": "/pipelines/smartseq2_single_sample/dockstore_SmartSeq2SingleSampleExample.json"
  },
  {
    "file_type": "SECONDARY_DESCRIPTOR",
    "path": "HISAT2.wdl"
  },
  {
    "file_type": "SECONDARY_DESCRIPTOR",
    "path": "Picard.wdl"
  },
  {
    "file_type": "SECONDARY_DESCRIPTOR",
    "path": "RSEM.wdl"
  }
]
```
Download each descriptor using GET /api/ga4gh/v2/tools/{id}/versions/{version_id}/{type}/descriptor/{relative_path} for both primary and secondary descriptors:
```
$ wget https://dockstore.org:8443/api/ga4gh/v2/tools/%23workflow%2Fgithub.com%2FHumanCellAtlas%2Fskylab/versions/master/PLAIN_WDL/descriptor//pipelines/smartseq2_single_sample/SmartSeq2SingleSample.wdl
$ wget https://dockstore.org:8443/api/ga4gh/v2/tools/%23workflow%2Fgithub.com%2FHumanCellAtlas%2Fskylab/versions/master/PLAIN_WDL/descriptor/HISAT2.wdl
$ wget https://dockstore.org:8443/api/ga4gh/v2/tools/%23workflow%2Fgithub.com%2FHumanCellAtlas%2Fskylab/versions/master/PLAIN_WDL/descriptor/Picard.wdl
$ wget https://dockstore.org:8443/api/ga4gh/v2/tools/%23workflow%2Fgithub.com%2FHumanCellAtlas%2Fskylab/versions/master/PLAIN_WDL/descriptor/RSEM.wdl
```
Download the test parameter file using GET /api/ga4gh/v2/tools/{id}/versions/{version_id}/{type}/tests:
```
$ wget https://dockstore.org:8443/api/ga4gh/v2/tools/%23workflow%2Fgithub.com%2FHumanCellAtlas%2Fskylab/versions/master/PLAIN_WDL/tests
```
To run the workflow locally using cromwell, you may first need to download all the files mentioned in the input parameter file and change it all the https path to local paths.  Then run the workflow using:
```
$ java -jar cromwell-*.jar run SmartSeq2SingleSample.wdl -i tests
```

## [TODO] Tagged releases
Pipelines are released in tagged releases to track changes in version of pipeline, and to track when a pipeline is blessed.

## [TODO] Automated testing options
## [TODO] Onboarding pipelines
