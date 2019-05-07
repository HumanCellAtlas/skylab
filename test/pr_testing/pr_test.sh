#!/usr/bin/env bash

# Find the git root
declare -r GIT_ROOT=$(git rev-parse --show-toplevel)

declare -r -a EXCLUDE_PIPELINES=(config pr_testing optimus smartseq2_single_sample)

function testDirShouldBeSkipped() {
  local -r testDir=${1}
  testDirName=$(basename ${testDir})
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
  for testDir in $(find ${GIT_ROOT}/test -type d -maxdepth 1 -mindepth 1); do
    testDirName=$(basename ${testDir})
    if $(testDirShouldBeSkipped ${testDir}); then
      continue
    fi
    echo

  done
}

findPipelinesToTest