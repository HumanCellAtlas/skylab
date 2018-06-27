#!/usr/bin/env python
# coding: utf-8
#
# EPY: stripped_notebook: {"metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"codemirror_mode": {"name": "ipython", "version": 3}, "file_extension": ".py", "mimetype": "text/x-python", "name": "python", "nbconvert_exporter": "python", "pygments_lexer": "ipython3", "version": "3.6.4"}}, "nbformat": 4, "nbformat_minor": 2}

# EPY: START code
# EPY: ESCAPE %matplotlib inline

import pysam
import itertools
import numpy as np
import scipy.sparse as sp
import matplotlib.pylab as plt
import pandas as pd
import natsort
import os
from collections import namedtuple
from collections import Counter
import distance

from typing import List, Set, Tuple
# EPY: END code

# EPY: START code
# (from https://samtools.github.io/hts-specs/SAMv1.pdf)
MASK_MULTI_SEG = 0x1
MASK_PROPER = 0x2
MASK_UNMAPPED = 0x4
MASK_NEXT_TMPL_UNMAPPED = 0x8
MASK_RC = 0x10
MASK_NEXT_RC = 0x20
MASK_FIRST_SEG = 0x40
MASK_LAST_SEG = 0x80
MASK_SA = 0x100
MASK_FAIL = 0x200
MASK_DUPL = 0x400
MASK_SUPL = 0x800
MASK_PRIMARY = 0x900
# EPY: END code

# EPY: START markdown
# ## Detailed analysis of discordant records
# EPY: END markdown

# EPY: START code
# `samtools sort -n merged.bam` 
op_bam_file = './optimus/merged.ns.bam'

# `samtools sort -n possorted_genome_bam.bam`
cr_bam_file = './cell_ranger/cellranger.ns.bam'

op_bam = pysam.AlignmentFile(op_bam_file)
cr_bam = pysam.AlignmentFile(cr_bam_file)
# EPY: END code

# EPY: START code
bam_record_attributes = namedtuple('bam_record_attributes', 'contig, start, end, cigar, flag, name')
bam_tags = namedtuple('bam_tags', 'CR, UR, CB')

def alignments_grouped_by_read_name_generator(bam_file: pysam.libcalignmentfile.AlignmentFile):
    """Iterates through a read-name-sorted BAM file and groups all alignments of a read.
    
    Returns:
        a tuple of read-name and a list of its alignments
    """
    bam_file.reset()
    for alignment in itertools.groupby(bam_file, key=lambda alignment: alignment.query_name):
        read_name = alignment[0]
        grouper = alignment[1]
        alignments = []
        try:
            while True:
                alignment = grouper.__next__()
                alignments.append(alignment)
        except StopIteration:
            pass
        yield read_name, alignments

def get_barcode_from_record(rec, barcode_tag_key='CR'):
    orig_barcode = None
    try:
        orig_barcode = rec.get_tag(barcode_tag_key)
    except:
        return None
    
    # record "-x" (for CellRanger)
    fixed_barcode = orig_barcode.split('-')[0]
    
    return fixed_barcode

def get_record_attributes(recs: List[pysam.AlignedSegment]) -> List[bam_record_attributes]:
    return [bam_record_attributes(
        contig=rec.reference_name,
        start=rec.reference_start,
        end=rec.reference_end,
        cigar=rec.cigarstring,
        flag=rec.flag,
        name=rec.query_name) for rec in recs]

def get_record_tags(recs: List[pysam.AlignedSegment]) -> List[bam_tags]:
    return [bam_tags(
        CR=get_barcode_from_record(rec, 'CR'),
        UR=get_barcode_from_record(rec, 'UR'),
        CB=get_barcode_from_record(rec, 'CB')) for rec in recs]

def get_sorted_attributes_tags_list(attr_list, tags_list):
    # sort by alignment
    sorted_index_attr_list = sorted(
        enumerate(attr_list),
        key=lambda entry: (entry[1].contig, entry[1].start, entry[1].end, entry[1].cigar))
    sort_index = [entry[0] for entry in sorted_index_attr_list]
    sorted_attr_list = [entry[1] for entry in sorted_index_attr_list]
    sorted_tags_list = [tags_list[j] for j in sort_index]
    return sorted_attr_list, sorted_tags_list
# EPY: END code

# EPY: START code
# instantiate generators for alignements grouped by read-name
op_grouped_records_gen = alignments_grouped_by_read_name_generator(op_bam)
cr_grouped_records_gen = alignments_grouped_by_read_name_generator(cr_bam)

# master list of parsed records (for further analysis)
op_master_rec_attr_list = list()
cr_master_rec_attr_list = list()
op_master_rec_tags_list = list()
cr_master_rec_tags_list = list()

# indices of discordant records
discordant_alignment_indices = list()
discordant_flag_indices = list()
discordant_CR_indices = list()
discordant_UR_indices = list()
discordant_CB_indices = list()

read_index = 0
for (op_read_name, op_recs), (cr_read_name, cr_recs) in zip(op_grouped_records_gen, cr_grouped_records_gen):
    assert op_read_name == cr_read_name, "Different read names"
    assert len(op_recs) == len(cr_recs), "Different number of alignment positions"
    
    op_rec_attr_list = get_record_attributes(op_recs)
    cr_rec_attr_list = get_record_attributes(cr_recs)
    op_tags_list = get_record_tags(op_recs)
    cr_tags_list = get_record_tags(cr_recs)
    
    # sort by alignment position
    sorted_op_rec_attr_list, sorted_op_tags_list = get_sorted_attributes_tags_list(
        op_rec_attr_list, op_tags_list)
    sorted_cr_rec_attr_list, sorted_cr_tags_list = get_sorted_attributes_tags_list(
        cr_rec_attr_list, cr_tags_list)
    
    # add to the master list
    op_master_rec_attr_list.append(sorted_op_rec_attr_list)
    cr_master_rec_attr_list.append(sorted_cr_rec_attr_list)
    op_master_rec_tags_list.append(sorted_op_tags_list)
    cr_master_rec_tags_list.append(sorted_cr_tags_list)
    
    # discordant records
    op_alignments = [(attr.contig, attr.start, attr.end, attr.cigar) for attr in sorted_op_rec_attr_list]
    cr_alignments = [(attr.contig, attr.start, attr.end, attr.cigar) for attr in sorted_cr_rec_attr_list]
    if op_alignments != cr_alignments:
        discordant_alignment_indices.append(read_index)

    op_flags = [attr.flag for attr in sorted_op_rec_attr_list]
    cr_flags = [attr.flag for attr in sorted_cr_rec_attr_list]
    if op_flags != cr_flags:
        discordant_flag_indices.append(read_index)

    op_CR = [tags.CR for tags in sorted_op_tags_list]
    cr_CR = [tags.CR for tags in sorted_cr_tags_list]
    if op_CR != cr_CR:
        discordant_CR_indices.append(read_index)

    op_UR = [tags.UR for tags in sorted_op_tags_list]
    cr_UR = [tags.UR for tags in sorted_cr_tags_list]
    if op_UR != cr_UR:
        discordant_UR_indices.append(read_index)

    op_CB = [tags.CB for tags in sorted_op_tags_list]
    cr_CB = [tags.CB for tags in sorted_cr_tags_list]
    if op_CB != cr_CB:
        discordant_CB_indices.append(read_index)

    # increment index
    read_index += 1
    
num_total_records = read_index
# EPY: END code

# EPY: START code
print(f'total number of reads: {num_total_records}')
print(f'number of reads with discordant alignments: {len(discordant_alignment_indices)}')
print(f'number of reads with discordant flags: {len(discordant_flag_indices)}, fraction = {len(discordant_flag_indices)/num_total_records:.4f}')
print(f'number of reads with discordant CR: {len(discordant_CR_indices)}, fraction = {len(discordant_CR_indices)/num_total_records:.4f}')
print(f'number of reads with discordant UR: {len(discordant_UR_indices)}, fraction = {len(discordant_UR_indices)/num_total_records:.4f}')
print(f'number of reads with discordant CB: {len(discordant_CB_indices)}, fraction = {len(discordant_CB_indices)/num_total_records:.4f}')
# EPY: END code

# EPY: START markdown
# ## Reads with discordant CB
# EPY: END markdown

# EPY: START code
op_missing_CB = 0
cr_missing_CB = 0

# list of tuples of d(op_CB, op_CR), d(cr_CB, cr_CR), d(op_CB, cr_CB)
hamming_distance_tuples = list()
implicated_CR_set = set()

for idx in discordant_CB_indices:
    op_CB_list = list(set(tags.CB for tags in op_master_rec_tags_list[idx]))
    cr_CB_list = list(set(tags.CB for tags in cr_master_rec_tags_list[idx]))
    assert len(op_CB_list) == 1
    assert len(cr_CB_list) == 1
    op_CB = op_CB_list[0]
    cr_CB = cr_CB_list[0]
    
    if op_CB is None:
        op_missing_CB += 1
        continue
        
    if cr_CB is None:
        cr_missing_CB += 1
        continue

    op_CR = op_master_rec_tags_list[idx][0].CR
    cr_CR = cr_master_rec_tags_list[idx][0].CR
    
    assert op_CR == cr_CR
    
    implicated_CR_set.add(op_CR)
    
    op_CB_CR_hamming_dist = distance.hamming(op_CB, op_CR)
    cr_CB_CR_hamming_dist = distance.hamming(cr_CB, cr_CR)    
    op_cr_CB_CB_hamming_dist = distance.hamming(op_CB, cr_CB)
    
    hamming_distance_tuples.append((op_CB_CR_hamming_dist, cr_CB_CR_hamming_dist, op_cr_CB_CB_hamming_dist))
    
op_CB_CR_hamming_dist_hist = Counter([entry[0] for entry in hamming_distance_tuples])
cr_CB_CR_hamming_dist_hist = Counter([entry[1] for entry in hamming_distance_tuples])
op_cr_CB_CB_hamming_dist_hist = Counter([entry[2] for entry in hamming_distance_tuples])

print(f'Optimus reads missing CB: {op_missing_CB}, fraction = {op_missing_CB/num_total_records:.4f}')
print(f'CellRanger reads missing CB: {cr_missing_CB}, fraction = {cr_missing_CB/num_total_records:.4f}')
print(f'Discordant CB between CellRanger and Optimus: {len(hamming_distance_tuples)}, fraction = {len(hamming_distance_tuples)/num_total_records:.4f}')
print(f'Distribution of (CB, CR) Hamming distance for Optimus reads: {op_CB_CR_hamming_dist_hist}')
print(f'Distribution of (CB, CR) Hamming distance for CellRanger reads: {cr_CB_CR_hamming_dist_hist}')
print(f'Distribution of (CB, CB) Hamming distance between Optimus and CellRanger reads: {op_cr_CB_CB_hamming_dist_hist}')
print(f'Number of unique implicated CRs: {len(implicated_CR_set)}')
# EPY: END code

# EPY: START markdown
# **Summary**
# 
# - All Optimus reads have CB tags
# - Some CellRanger reads do not have CB tags
# - Both Optimis and CellRanger only correct barcodes within 1HD from the whitelist
# - Optimus and CellRanger correct barcode errors in _different ways_
# EPY: END markdown

# EPY: START markdown
# ## Reads with discordant flags
# EPY: END markdown

# EPY: START code
class FlagSummaryStatistics:
    def __init__(self):
        self.num_multi_seg = 0
        self.num_proper = 0
        self.num_umapped = 0
        self.num_next_templ_unmapped = 0
        self.num_rc = 0
        self.num_next_rc = 0
        self.num_first_seg = 0
        self.num_last_seg = 0
        self.num_sa = 0
        self.num_fail = 0
        self.num_dupl = 0
        self.num_supl = 0
    
    def _add_flag_to_summary(self, flag):
        self.num_multi_seg += (flag & MASK_MULTI_SEG) >> 0
        self.num_proper += (flag & MASK_PROPER) >> 1
        self.num_umapped += (flag & MASK_UNMAPPED) >> 2
        self.num_next_templ_unmapped += (flag & MASK_NEXT_TMPL_UNMAPPED) >> 3
        self.num_rc += (flag & MASK_RC) >> 4
        self.num_next_rc += (flag & MASK_NEXT_RC) >> 5
        self.num_first_seg += (flag & MASK_FIRST_SEG) >> 6
        self.num_last_seg += (flag & MASK_LAST_SEG) >> 7
        self.num_sa += (flag & MASK_SA) >> 8
        self.num_fail += (flag & MASK_FAIL) >> 9
        self.num_dupl += (flag & MASK_DUPL) >> 10
        self.num_supl += (flag & MASK_SUPL) >> 11
        
    @staticmethod
    def _combine_flags(flags, op='or', primary_only=False):
        if primary_only:
            try:
                filtered_flags = [flag for flag in flags if flag & MASK_PRIMARY == 0]
            except:
                raise Exception("No primary alignment is available")
        else:
            filtered_flags = flags
            
        combined_flag = 0
        if op == 'or':
            for flag in filtered_flags:
                combined_flag |= flag
        elif op == 'xor':
            for flag in filtered_flags:
                combined_flag ^= flag
        elif op == 'and':
            for flag in filtered_flags:
                combined_flag &= flag
        else:
            raise Exception(f'Unknown bitwise operation {op}')

        return combined_flag
    
    def process_single(self, flags, op='or', primary_only=False):
        """Takes an iterable of binary flags, combines them according to a given bitwise operation, and updates
        the summary statistics."""
        combined_flag = self._combine_flags(flags, op, primary_only)
        self._add_flag_to_summary(combined_flag)

    def process_difference(self, first_flags, second_flags, joint_op='xor', primary_only=False):
        first_combined_flag = self._combine_flags(first_flags, op='or', primary_only=primary_only)
        second_combined_flag = self._combine_flags(second_flags, op='or', primary_only=primary_only)
        if joint_op == 'xor':
            self._add_flag_to_summary(first_combined_flag ^ second_combined_flag)
        elif joint_op == 'and':
            self._add_flag_to_summary(first_combined_flag & second_combined_flag)
        else:
            raise Exception(f'Unknown joint operation {joint_op}')
            
    def print_summary(self):
        print(f'num_multi_seg: {self.num_multi_seg}')
        print(f'num_proper: {self.num_proper}')
        print(f'num_umapped: {self.num_umapped}')
        print(f'num_next_templ_unmapped: {self.num_next_templ_unmapped}')
        print(f'num_rc: {self.num_rc}')
        print(f'num_next_rc: {self.num_next_rc}')
        print(f'num_first_seg: {self.num_first_seg}')
        print(f'num_last_seg: {self.num_last_seg}')
        print(f'num_sa: {self.num_sa}')
        print(f'num_fail: {self.num_fail}')
        print(f'num_dupl: {self.num_dupl}')
        print(f'num_supl: {self.num_supl}')
# EPY: END code

# EPY: START code
op_multiple_primary = 0
cr_multiple_primary = 0

for idx in discordant_flag_indices:
    op_flags = [attr.flag for attr in op_master_rec_attr_list[idx]]
    cr_flags = [attr.flag for attr in cr_master_rec_attr_list[idx]]
    
    # reads with more than one primary alignment (!)
    # according to SAM specs, each read should have only one primary alignment
    # Note: the issue stems from the synthetic data -- will follow up; for now, let us
    # ignore such reads
    op_recs_multiple_primary = sum(flag & MASK_PRIMARY == 0 for flag in op_flags) > 1
    cr_recs_multiple_primary = sum(flag & MASK_PRIMARY == 0 for flag in cr_flags) > 1
    op_multiple_primary += op_recs_multiple_primary
    cr_multiple_primary += cr_recs_multiple_primary
# EPY: END code

# EPY: START code
op_flag_summary_stats = FlagSummaryStatistics()
cr_flag_summary_stats = FlagSummaryStatistics()

for idx in discordant_flag_indices:
    op_flags = [attr.flag for attr in op_master_rec_attr_list[idx]]
    cr_flags = [attr.flag for attr in cr_master_rec_attr_list[idx]]

    op_flag_summary_stats.process_single(op_flags, primary_only=False)
    cr_flag_summary_stats.process_single(cr_flags, primary_only=False)
# EPY: END code

# EPY: START code
op_primary_only_flag_summary_stats = FlagSummaryStatistics()
cr_primary_only_flag_summary_stats = FlagSummaryStatistics()

for idx in discordant_flag_indices:
    op_flags = [attr.flag for attr in op_master_rec_attr_list[idx]]
    cr_flags = [attr.flag for attr in cr_master_rec_attr_list[idx]]

    op_primary_only_flag_summary_stats.process_single(op_flags, primary_only=True)
    cr_primary_only_flag_summary_stats.process_single(cr_flags, primary_only=True)
# EPY: END code

# EPY: START code
primary_only_flag_xor_summary_stats = FlagSummaryStatistics()
primary_only_flag_and_summary_stats = FlagSummaryStatistics()

for idx in discordant_flag_indices:
    op_flags = [attr.flag for attr in op_master_rec_attr_list[idx]]
    cr_flags = [attr.flag for attr in cr_master_rec_attr_list[idx]]

    primary_only_flag_xor_summary_stats.process_difference(
        op_flags, cr_flags, primary_only=True, joint_op='xor')

    primary_only_flag_and_summary_stats.process_difference(
        op_flags, cr_flags, primary_only=True, joint_op='and')
# EPY: END code

# EPY: START code
print(f'Number of Optimus reads with multiple primary alignments: {op_multiple_primary}, fraction = {op_multiple_primary/num_total_records:.4f}')
print(f'Number of CellRanger reads with multiple primary alignments: {cr_multiple_primary}, fraction = {cr_multiple_primary/num_total_records:.4f}')
print()

print('Optimus (primary only):')
op_primary_only_flag_summary_stats.print_summary()
print()

print('CellRanger (primary only):')
cr_primary_only_flag_summary_stats.print_summary()
print()

print('joint XOR (primary only):')
primary_only_flag_xor_summary_stats.print_summary()
print()

print('joint AND (primary only):')
primary_only_flag_and_summary_stats.print_summary()
print()

print('Optimus (all):')
op_flag_summary_stats.print_summary()
print()

print('CellRanger (all):')
cr_flag_summary_stats.print_summary()
print()
# EPY: END code

# EPY: START markdown
# **Summary**
# 
# - There are about ~ 70k identical raw reads in the synthetic data and they need to be removed
# - For primary alignments with discordant flags, 47850 reads are marked as duplicate by both Optimus and CellRanger
#   whereas 446777 reads (about half of all reads) are marked duplicate differently
# - CellRanger marks more primary alignments as duplicate (356877 vs 185600), however, this does not necessary
#   imply more molecules.
# - It looks like `sctools` does not use the duplicate flag at all; `MarkDuplicates` only _marks_ duplicates. It does not remove them.
# EPY: END markdown

# EPY: START code
num_primary = 0
for record in pysam.AlignmentFile("./optimus/cell-sorted.bam"):
    if (record.flag & MASK_PRIMARY) == 0 and record.mapping_quality == 255:
        num_primary += 1
# EPY: END code

# EPY: START code
num_primary
# EPY: END code

# EPY: START code
num_primary = 0
for record in pysam.AlignmentFile("./cell_ranger/possorted_genome_bam.bam"):
    if (record.flag & MASK_PRIMARY) == 0 and record.mapping_quality == 255:
        num_primary += 1
# EPY: END code

# EPY: START code
num_primary
# EPY: END code

# EPY: START markdown
# **Summary**
# 
# -- CellRanger definitely fiddles with flags in non-standard ways
# EPY: END markdown

# EPY: START markdown
# ## A brute-force expression matrix calculator
# 
# The idea is to mimick CellRanger.
# EPY: END markdown

# EPY: START code
from sctools import gtf
from typing import List, Dict, Tuple, Set
from scipy import sparse as sp
import operator
# EPY: END code

# EPY: START code
sorted_n_GE_UB_CB_bam_file = './optimus/sorted_n_GE_UB_CB.bam'
annotation_file = './references/cellranger.gtf'
gene_name_tag = 'GE'
# EPY: END code

# EPY: START code
gene_name_to_index: Dict[str, int] = {}
gtf_reader = gtf.Reader(annotation_file)

# map the gene from reach record to an index in the sparse matrix
for gene_index, record in enumerate(gtf_reader.filter(retain_types=['gene'])):
    gene_name = record.get_attribute('gene_name')
    if gene_name is None:
        raise ValueError(
            'malformed GTF file detected. Record is of type gene but does not have a '
            '"gene_name" field: %s' % repr(record))
    gene_name_to_index[gene_name] = gene_index

n_genes = len(gene_name_to_index)
# EPY: END code

# EPY: START code
def alignments_grouped_by_query_name_generator(bam_path: pysam.libcalignmentfile.AlignmentFile,
                                               cell_barcode_tag: str,
                                               molecule_barcode_tag: str):
    """Iterates through a query-name-sorted BAM file and groups all alignments of a read.
    
    Returns:
        a tuple of read-name and a list of its alignments
    """
    with pysam.AlignmentFile(bam_path, 'rb') as bam_file:
        for alignment in itertools.groupby(bam_file, key=lambda alignment: alignment.query_name):
            query_name = alignment[0]
            grouper = alignment[1]
            alignments = []
            try:
                while True:
                    alignment = grouper.__next__()
                    alignments.append(alignment)
            except StopIteration:
                pass

            cell_barcode = None
            try:
                cell_barcode = alignments[0].get_tag(cell_barcode_tag)
            except:
                pass

            molecule_barcode = None
            try:
                molecule_barcode = alignments[0].get_tag(molecule_barcode_tag)
            except:
                pass

            yield query_name, cell_barcode, molecule_barcode, alignments
# EPY: END code

# EPY: START code
cell_barcode_tag = 'CB'
molecule_barcode_tag = 'UB'
gene_name_tag = 'GE'

# keep track of the observed triples (cell_barcode, molecule_barcode, gene_name) in a hash set
observed_cell_molecule_gene_set: Set[Tuple[str, str, str]] = set()
    
# COO sparse matrix entries
data: List[int] = []
cell_indices: List[int] = []
gene_indices: List[int] = []
    
# track which cells we've seen, and what the current cell number is
n_cells = 0
cell_barcode_to_index: Dict[str, int] = {}

grouped_records_generator = alignments_grouped_by_query_name_generator(
    sorted_n_GE_UB_CB_bam_file, cell_barcode_tag, molecule_barcode_tag)

for query_name, cell_barcode, molecule_barcode, alignments in grouped_records_generator:
    
    if cell_barcode is None or molecule_barcode is None: # only keep queries w/ well-formed UMIs
        continue
    
    if len(alignments) == 1:
        primary_alignment = alignments[0]
        if primary_alignment.has_tag(gene_name_tag):
            gene_name = primary_alignment.get_tag(gene_name_tag)
        else:
            continue # drop query
    else: # multi-map
        implicated_gene_names: Set[str] = set()
        for alignment in alignments:
            if alignment.has_tag(gene_name_tag):
                implicated_gene_names.add(alignment.get_tag(gene_name_tag))
        if len(implicated_gene_names) == 1: # only one gene
            gene_name = implicated_gene_names.__iter__().__next__()
        else:
            continue # drop query
    
    if (cell_barcode, molecule_barcode, gene_name) in observed_cell_molecule_gene_set:
        continue # optical/PCR duplicate -> drop query
    else:
        observed_cell_molecule_gene_set.add((cell_barcode, molecule_barcode, gene_name))
    
    # find the indices that this molecule should correspond to
    gene_index = gene_name_to_index[gene_name]

    # if we've seen this cell before, get its index, else set it
    try:
        cell_index = cell_barcode_to_index[cell_barcode]
    except KeyError:
        cell_index = n_cells
        cell_barcode_to_index[cell_barcode] = n_cells
        n_cells += 1
        
    # record the molecule data
    data.append(1)  # one count of this molecule
    cell_indices.append(cell_index)
    gene_indices.append(gene_index)
    
# convert into coo_matrix
coordinate_matrix = sp.coo_matrix((data, (cell_indices, gene_indices)),
    shape=(n_cells, n_genes), dtype=np.uint32)

# convert into csr matrix and return
col_index = np.asarray([k for k, v in sorted(gene_name_to_index.items(), key=operator.itemgetter(1))])
row_index = np.asarray([k for k, v in sorted(cell_barcode_to_index.items(), key=operator.itemgetter(1))])
# EPY: END code

# EPY: START code
print(f'optimus total counts: {np.sum(coordinate_matrix)}')
# EPY: END code

# EPY: START code
import sctools
import collections
import tables
import scipy.sparse as sp_sparse

GeneBCMatrix = collections.namedtuple('GeneBCMatrix', ['gene_ids', 'gene_names', 'barcodes', 'matrix'])
 
def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print("That genome does not exist in this file.")
            return None
        gene_ids = getattr(group, 'genes').read()
        gene_names = getattr(group, 'gene_names').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return GeneBCMatrix(gene_ids, gene_names, barcodes, matrix)
 
filtered_matrix_h5 = "cell_ranger/raw_gene_bc_matrices_h5.h5"
genome = "GRCh38"
cr_cm = get_matrix_from_h5(filtered_matrix_h5, genome)
# EPY: END code

# EPY: START code
print(f'cell ranger total counts: {cr_cm.matrix.sum()}')
# EPY: END code

# EPY: START markdown
# **Summary**
# 
# - Close enough -- 95% concordance.
# - I do not quite understand the difference; some thoughts:
#   - Duplicate marking is immaterial for 10x since we only count unique (CB, GE, UB) tuples.
#   - Multiple primary alignment is not an issue (unles CellRanger throws away reads with multiple primary alignments?)
#   - Filtering by MQ could be an issue (Optimus does not filter, CellRanger does, but it manipulates the STAR bam out such that it is difficult to interpret their downstream filter MQ == 255)
#   - Another source of discrepancy is the difference in barcode error correction
#   - A thorough reading of CellRanger code (WIP) is helpful.
# EPY: END markdown
