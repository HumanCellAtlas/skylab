#!/usr/bin/env python3

import click

@click.command()
@click.option('--in-fastq-gz',help='gz compressed input fastq file')
@click.option('--out-fastq-gz',help='gz compressed output fastq file')
@click.option('--keep-reads-gz',help='gz compressed file with read names to keep')
@click.option('--verbose',help='verbose',default=False,is_flag=True,flag_value=True)
def filterbyreadname(in_fastq_gz, out_fastq_gz, keep_reads_gz, verbose):
    import gzip
    import sys
    from Bio import SeqIO
    ## Put reads that are in the keep list in the output file
    keep_read_set = set(line.decode('ascii').rstrip() for line in gzip.open(keep_reads_gz,'rb'))
    if (verbose):
        print('Done loading keep read list', file=sys.stderr);
    ## Loop over input and filter reads
    with gzip.open(out_fastq_gz, 'wt') as output_file:
        with gzip.open(in_fastq_gz, 'rt') as input_file:
            counter=0;
            for rec in SeqIO.parse(input_file, 'fastq'):
                counter+=1;
                if(verbose and counter % 1e5 == 0):
                    print('Processed {} reads'.format(counter));
                if(rec.id in keep_read_set):
                    SeqIO.write(rec, output_file, 'fastq');
    ## Print completed message
    if(verbose):
        print('Completed',file=sys.stderr)

if __name__ == '__main__':
    filterbyreadname()
