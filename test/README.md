**Scientific and PR testing infrastructure in Skylab.**

Both tests share the same infrastructure, and follow the portability spec; the WDLs in this PR represent a cromwell test case for the [portability test](https://docs.google.com/document/d/1ghLoHMbKOPsndA1WgdSAHm5X82p86ryLBiAt1hz6HuI/edit); they are just missing a docker environment to run the WDLs in. 

The PR test is designed to answer the question: "did my changes result in any change to the outputs of my pipeline?"

The scientific test is designed to answer the question: "Is the pipeline performing within an expected range of scientific outputs?"

**Testing Logic & Files**
All files added in this PR are contained within the test directory. 
- `tests/trigger_test.sh` triggers tests to run within the `humancellatlas/cromwell-tools` docker. This is the test called by Jenkins, and it takes parameters: inputs, wdl, and dependencies. It is designed to trigger tests in the form of WDL workflows. 
- `tests/test_cromwell_workflow.sh` is called by `tests/trigger_test.sh` within the cromwell-tools docker. It calls the provided wdl, waits for it to finish, and checks it's success state. 
- The scientific and PR tests both contain an "infrastructure wdl" which follows the form `test_${PIPELINE_FOLDER_NAME}_PR` and is composed of two parts:
  - The pipeline workflow, which is imported and run to completion with call caching, which allows the test to run quickly when nothing has changed that would modify the pipeline. 
  - The checker workflow, which, given the outputs of the just-run pipeline, tests them either against the expected outputs (PR) or a range of acceptable outputs (scientific) 

The tests are designed so that any failure results in `exit 1` within the workflow, after the reason is logged to stderr. If a workflow fails due to a result not matching expectations, the result will be logged in the stderr produced by the checker task. If a workflow fails due to not running properly, the user will need to go find the part of the workflow that failed. The jenkins job and cromwell backend will both log the workflow ID so that it is easy to go back and find the failed workflow. 

**Updating workflows that result in scientific changes:**
In the event that a PR produces a desirable scientific change, the PR test will fail. In order to ensure that future PRs test against the correct "expected values", the test inputs must be updated to reflect the success state of the incoming pipeline. For optimus, this means updating the md5 hashes corresponding to the new output files, which are found in: 
`test/optimus/pr/test_inputs.json` (PR)
`test/optimus/scientific/test_inputs.json` (scientific)

More complex tests may require different inputs, but the expected values should always be parameterized in the `test_inputs.json` files to facilitate easy updating. 

**Jenkins Testing details**

The jenkins job executes as follows: 
```bash
VAULT_TOKEN=$(cat /etc/vault-token-dsde)
PIPELINE_FOLDER_NAME=<optimus|smartseq2_single_sample>
INPUTS_JSON="https://raw.githubusercontent.com/HumanCellAtlas/skylab/master/test/${PIPELINE_FOLDER_NAME}/pr/test_inputs.json"
WDL_FILE="/working/test/${PIPELINE_FOLDER_NAME}/pr/Test${PIPELINE_FOLDER_NAME}_PR.wdl"
DEPENDENCIES_JSON="/working/test/${PIPELINE_FOLDER_NAME}/pr/dependencies.json"
bash ./test/trigger_test.sh $VAULT_TOKEN $INPUTS_JSON $WDL_FILE $DEPENDENCIES_JSON
```

Thus, `trigger_test.sh` is run from the current branch, which is cloned as the working directory
of the jenkins instance. 
`trigger_test.sh` then runs a docker that maps the current working directory to `/working`.
Thus, the `WDL_FILE`, and `DEPENDENCIES_JSON` (and any files inside the dependencies with local paths) are interpreted from within the docker as referencing the relevant files on the current branch. 