#!/usr/bin/env python

## Author: Nick Barkas <nbarkas@broadinstitute.org>
## Date: May 9th, 2019
## Description: Pre-process fastq files by attaching barcodes to read names (generally a bad idea)

import argparse
from Bio import SeqIO
import time
import gzip

parser = argparse.ArgumentParser(description="Extract barcodes from index reads into read names for scATAC-seq BICCN data V1")
parser.add_argument('--index1', dest='file_i1', help='index1 FASTQ file')
parser.add_argument('--index2',dest='file_i2',help='index2 FASTQ file')
parser.add_argument('--read1', dest='file_r1',help='read1 FASTQ file')
parser.add_argument('--read2',dest='file_r2',help='read2 FASTQ files')
parser.add_argument('--out-read1', dest='out_file_r1', help='output file read1')
parser.add_argument('--out-read2', dest='out_file_r2', help='output file read2')

args = parser.parse_args()

# Open the files
fh_i1 = gzip.open(args.file_i1, 'rU')
fh_i2 = gzip.open(args.file_i2, 'rU')
fh_r1 = gzip.open(args.file_r1, 'rU')
fh_r2 = gzip.open(args.file_r2, 'rU')
# output files
fh_out_r1 = gzip.open(args.out_file_r1, 'w')
fh_out_r2 = gzip.open(args.out_file_r2, 'w')

# Get iterators for the files
iter_i1 = SeqIO.parse(fh_i1, 'fastq')
iter_i2 = SeqIO.parse(fh_i2, 'fastq')
iter_r1 = SeqIO.parse(fh_r1, 'fastq')
iter_r2 = SeqIO.parse(fh_r2, 'fastq')

start_time = time.time()

nEntries = 0
readComplete = False;
while not readComplete:
    if (nEntries % 100000 == 0):
        print("Processed %i records..." % nEntries)
    try:
        i1_sreq = iter_i1.next();
        i2_sreq = iter_i2.next();
        r1_sreq = iter_r1.next();
        r2_sreq = iter_r2.next();

        # Extract CB
        i1a = i1_sreq.seq[1:8]
        i1b = i1_sreq.seq[36:43]
        i2a = i2_sreq.seq[1:8]
        i2b = i2_sreq.seq[30:37]
        cb = str(i1a + i1b + i2a + i2b)

        # Append CB to read names
        r1_newname = cb + ':' + r1_sreq.name
        r2_newname = cb + ':' + r2_sreq.name

        # Create new read entries with the amended names
        r1_sreq_new = SeqIO.SeqRecord(r1_sreq.seq, r1_newname, r1_newname, r1_sreq.description, r1_sreq.dbxrefs,
                                      r1_sreq.features, r1_sreq.annotations, r1_sreq.letter_annotations)
        r2_sreq_new = SeqIO.SeqRecord(r2_sreq.seq, r2_newname, r2_newname, r2_sreq.description, r2_sreq.dbxrefs,
                                      r2_sreq.features, r2_sreq.annotations, r2_sreq.letter_annotations)

        # Write output
        SeqIO.write(r1_sreq_new, fh_out_r1, 'fastq')
        SeqIO.write(r2_sreq_new, fh_out_r2, 'fastq')

        nEntries += 1

    except StopIteration:
        print('done reading files')
        readComplete = True;

        ## TODO: Check if all the other iterators do not 'have_next()' if they do we need to issue a user warning


end_time = time.time()
print('Processed %i entries in %d seconds' % (nEntries, end_time - start_time))

# Close the input fh
fh_i1.close()
fh_i2.close()
fh_r1.close()
fh_r2.close()

# Close the output filehandles
fh_out_r1.close()
fh_out_r2.close()



