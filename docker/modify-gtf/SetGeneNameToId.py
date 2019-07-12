#!/usr/bin/env python

import argparse
import re

parser = argparse.ArgumentParser(description="Set the gene_name within a gtf to be equivalent to the values within gene_id.")
parser.add_argument('--in-gtf-file', dest='ingtf', help='input gtf file')
parser.add_argument('--out-gtf-file', dest='outgtf', help='output gtf file')

args = parser.parse_args()


def setGeneNameToId(in_gtf, out_gtf, verbose=True):
    with open(in_gtf, 'r') as fpin, open(out_gtf, 'w') as fpout:
        for line in fpin:
            stripped_line = line.strip()
            gene_id = re.search(r'gene_id ([^;]*);', stripped_line)
            gene_name = re.search(r'gene_name ([^;]*);', stripped_line)
            if gene_id and gene_name:
                modified_line = re.sub(r'gene_name ([^;]*);', 'gene_name ' + gene_id.group(1) + ";", stripped_line)
                fpout.write(modified_line + '\n')
            else:
                fpout.write(stripped_line + '\n')

setGeneNameToId(args.ingtf, args.outgtf)
