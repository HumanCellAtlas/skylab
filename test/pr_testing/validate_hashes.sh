#!/usr/bin/env bash

# traversal is for other people
declare -r GIT_ROOT=$(git rev-parse --show-toplevel)
declare -r HASHES_FILE=${GIT_ROOT}/test/test_a/pr/TestA.hashes.json
declare -a KEYS=$(jq -r 'keys' ${HASHES_FILE})
declare -r HASH="0cc175b9c0f1b6a831c399e269772661"
for KEY in $(echo "${KEYS}" | jq -r '.[]')
do
  # need `jq -r` to get results without surrounding double-quotes
  if [[ $HASH = $(jq -r '.["'$KEY'"]' $HASHES_FILE) ]]
  then
    echo "yay"
  else
    echo "HASH: $HASH"
    echo "HASHES_FILE_VALUE: $(jq -r '.["'$KEY'"]' $HASHES_FILE)"
  fi
done