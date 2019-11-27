#!/usr/bin/env bash

# Usage function
function usage() {
  echo "Usage: $0 [-m] -i input_directory -o output_directory"
}

# Defaults

# Command to issue on each input zarr file,
# default to cp, but can be set to mv to conserve
# disk space and for faster operation
execCommand='cp'

# Read the CLI arguments
while getopts ":i:o:hm" opt; do
  case "${opt}" in
    i)
      inputDir=$OPTARG
      ;;
    o)
      outputDir=$OPTARG
      ;;
    h)
      usage
      exit 0;
      ;;
    m)
      execCommand='mv'
      ;;
    *)
      usage
      exit 1;
      ;;
    esac
done
shift $((OPTIND -1))

# Check the inputs
if [ -z "$inputDir" ]; then
  echo "Error: input directory (-i) not specified"
  exit 1;
fi
if [ -z "$outputDir" ]; then
  echo "Error: output directory (-o) not specified"
  exit 1;
fi
if [ ! -d "$inputDir" ]; then
  echo "Error the input directory path is not a directory";
  exit 1;
fi
if [ ! -d "$outputDir" ]; then
  echo "Error; the output directory path is not a directory";
  exit;
fi

# Do the conversion
while IFS= read -r -d '' f
do
  bname1=$(basename "${f}")
  rel_final_path=$(echo "${bname1}" | tr '!' '/')
  rel_final_dirname=$(dirname "${rel_final_path}")
  rel_final_bname=$(basename "${rel_final_path}")
  mkdir -p "${outputDir}/${rel_final_dirname}"
  ${execCommand} "${f}" "${outputDir}/${rel_final_dirname}/${rel_final_bname}"
done <   <(find "${inputDir}" -type f -print0 )

