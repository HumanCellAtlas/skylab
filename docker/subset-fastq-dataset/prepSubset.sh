#!/bin/bash --

## Description: Given fastqs and their alignments (in bam format)
## filter the input fastqs according to their alignment position

set -eo pipefail

show_help() {
    echo "Usage: $0 [arguments]"
    echo ""
    echo " -a      google bucket location with the output of star-Align"
    echo " -f      google bucket with fastq input"
    echo " -r      region of alignments to keep"
    echo " -h      print this helpful message"
    echo " -o      output prefix"
    echo ""
}

## Process command line arguments
OPTIND=1

starAlignOutBucket=""
inFastqBucket=""
samtoolsRegion=""
outputPrefix=""

while getopts "ha:f:r:o:" opt; do
    case "$opt" in
	h)
	    show_help
	    exit 0
	    ;;
	a)
	    starAlignOutBucket=$OPTARG
	    ;;
	f)
	    inFastqBucket=$OPTARG
	    ;;
	r)
	    samtoolsRegion=$OPTARG
	    ;;
	o)
	    outputPrefix=$OPTARG
	    ;;
    esac
done

shift $((OPTIND-1))

## DEBUG check args
# echo starAlignOutBucket $starAlignOutBucket
# echo inFastqBucket $inFastqBucket
# echo samtoolsRegion $samtoolsRegion

## Input validation
if [ -z "$starAlignOutBucket" ]; then
    echo "Star align bucket location not provided"
    exit 1;
fi

if [ -z "$inFastqBucket" ]; then
    echo "Fastq bucket location not provided"
    exit 1;
fi

if [ -z "samtoolsRegion" ]; then
    echo "Genomic region to select not provided"
    exit 1;
fi

if [ -z "$outputPrefix" ]; then
    echo "Output prefix is empty"
    exit 1;
fi

###########################
## Internal params
keepreadsgz=keepreads.txt.gz
samtoolsSortCores=4
samtoolsMemPerCore=5G

## Get the aligned files
echo Downloading alignments...
mkdir starAlignOut
gsutil -m cp -r $starAlignOutBucket starAlignOut/

## Get the fastqs
echo Downloading fastqs...
mkdir fastqs
gsutil -m cp $inFastqBucket/* fastqs/

## Merge the bams into one file
echo Concatenating alignments...
samtools cat `find ./starAlignOut/ -name '*.bam' -printf '%p '` > concat.bam

## Sort and index the concatenated file
echo Sorting merged alignments...
samtools sort -@ ${samtoolsSortCores} -m ${samtoolsMemPerCore} concat.bam > concat.sorted.bam
samtools index concat.sorted.bam

## Extract read names
echo Extracting read names for $samtoolsRegion
samtools view concat.sorted.bam $samtoolsRegionn | cut -f 1 | pigz > $keepreadsgz

## Filter all the input fastqs in parallel
echo Filtering fastqs...
for infastqgz in `ls fastqs/*.fastq.gz`
do
    echo Processing ${infastqgz}
    outfastqgz=${outputPrefix}/`basename ${infastqgz}`
    filterFastqByReadName.py \
             --in-fastq-gz $infastqgz \
             --out-fastq-gz $outfastqgz \
             --keep-reads-gz $keepreadsgz \
             --verbose &
    sleep 10
done
wait;
