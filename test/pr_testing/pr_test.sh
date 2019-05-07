#!/usr/bin/env bash

# Find the git root
declare -r GIT_ROOT=$(git rev-parse --show-toplevel)

declare -r -a EXCLUDE_PIPELINES=(config pr_testing optimus smartseq2_single_sample)

function getPipelineDir() {
  local -r dirName=${1}
  echo ${GIT_ROOT}/pipelines/${dirName}
}

function getTestDir() {
  local -r dirName=${1}
  echo ${GIT_ROOT}/test/${dirName}
}

function testDirShouldBeSkipped() {
  local -r testDirName=${1}
  local -r testDir=$(getTestDir ${testDirName})
  if [[ " ${EXCLUDE_PIPELINES[*]} " == *"${testDirName}"* ]]; then
    >&2 echo "${testDir} is excluded from PR testing"
    echo true
  elif ! find ${testDir} -type d -name pr | grep '.*' > /dev/null; then
    >&2 echo "${testDir} does not have any PR tests"
    echo true
  else
    echo false
  fi
}

function findPipelinesToTest() {
  for testDirName in $(find ${GIT_ROOT}/test -type d -maxdepth 1 -mindepth 1 -exec basename {} \;); do

    # if this directory shouldn't be tested, skip it
    if $(testDirShouldBeSkipped ${testDirName}); then
      continue
    fi

    wdlDir=$(getPipelineDir ${testDirName})
    # there should only be one WDL in the pipelines dir



    wdlName=$(find ${wdlDir} -type f -name '*.wdl' -maxdepth 1 | head -n 1)

    >&2 echo ${wdlName} has a PR test

    echo ${testDirName}
  done
}

dirsToTest=($(findPipelinesToTest 2> /dev/null))

for dir in ${dirsToTest[@]}; do
  echo "${dir} has pr tests"
done