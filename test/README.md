**Running Tests Manually**
1) Set the environment variable BROAD_CROMWELL_KEY to have the contents of the cromwell credentials key JSON file.
2) cd to the top level directory of the repository
3) run the following command ```./test/trigger_test.sh [name of workflow without trailing /pr/]```

Example:
```
export BROAD_CROMWELL_KEY=`cat ~/identities/credentials.json`
cd ~/skylab/
./test/trigger_test.sh optimus
```

**PR testing infrastructure in Skylab.**
The tests follow the portability spec; the WDLs in this PR represent a cromwell test case for the [portability test](https://docs.google.com/document/d/1ghLoHMbKOPsndA1WgdSAHm5X82p86ryLBiAt1hz6HuI/edit); they are just missing a docker environment to run the WDLs in. 

The PR test is designed to answer the question: "did my changes result in any change to the outputs of my pipeline?"

**Testing Logic & Files**
All files related to testing are in the test/ repository
- `tests/trigger_test.sh` starts a `quay.io/broadinstitute/cromwell-tools` docker and launches tests on caas cromwell
- `tests/test_cromwell_workflow.sh` is called by `tests/trigger_test.sh` within the cromwell-tools docker. It calls the provided wdl, waits for it to finish, and checks it's success state. 
- The scientific and PR tests both contain an "infrastructure wdl" which follows the form `test_${PIPELINE_FOLDER_NAME}_PR.wdl` and is composed of two parts:

  - The pipeline workflow, which is imported and run to completion with call caching, which allows the test to run quickly when nothing has changed that would modify the pipeline. 
  - The checker workflow, which, given the outputs of the just-run pipeline, tests them either against the expected outputs (PR) or a range of acceptable outputs (scientific) 

The tests are designed so that any failure results in `exit 1` within the workflow, after the reason is logged to stderr. If a workflow fails due to a result not matching expectations, the result will be logged in the stderr produced by the checker task.

If a workflow fails due to not running properly, the user will need to go find the part of the workflow that failed. The jenkins job and cromwell backend will both log the workflow ID so that it is easy to go back and find the failed workflow.

**Updating workflows that result in scientific changes:**
In the event that a PR produces a desirable scientific change, the PR test will fail. In order to ensure that future PRs test against the correct "expected values", the test inputs must be updated to reflect the success state of the incoming pipeline. For optimus, this means updating the md5 hashes corresponding to the new output files, which are found in: 
`test/optimus/pr/test_inputs.json` (PR)

More complex tests may require different inputs, but the expected values should always be parameterized in the `test_inputs.json` files to facilitate easy updating. 
