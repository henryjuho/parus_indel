#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
import sys
import os
import gzip
sys.path.insert(0, os.getenv('HOME') + '/parus_indel/summary_analyses')
sys.path.insert(0, os.getenv('HOME') + '/parus_indel/anavar_analyses')
from indel_lengths import indel_type
from window_sel_vs_neu_anavar import window_call_sites


def missing_genotypes(vcf_line):

    """
    flags vcf lines with genotypes missing
    :param vcf_line: pysam.VariantRecord
    :return: bool
    """

    genotypes = [vcf_line.samples[x]['GT'] for x in vcf_line.samples.keys()]
    geno_list = sum([list(x) for x in genotypes], [])

    if None in geno_list:
        return True
    else:
        return False


def print_dict(len_dict, chromo, region_id, n_callable, include_header=False):

    """
    prints length dict
    :param len_dict: dict
    :param chromo: str
    :param region_id: str
    :param n_callable: int
    :param include_header: bool
    :return: None
    """

    indivs = sorted(len_dict.keys())
    indivs.remove('n_indel')
    indivs.remove('total_event_length')

    # prints header
    if include_header:
        header = ['chr'] + indivs + ['total', 'n_indel', 'total_event_length', 'indel_type', 'region', 'callable']
        print(*header, sep=',')

    # prints data
    for x in ['ins', 'del']:
        bp_list = [len_dict[z][x] for z in indivs]
        data_line = [chromo] + bp_list + [sum(bp_list), len_dict['n_indel'][x], len_dict['total_event_length'][x],
                                          x, region_id, n_callable]
        print(*data_line, sep=',')


def sum_loss_gain_dict(dict1, dict2):

    """
    adds two sequence length dictionaries
    :param dict1: dict
    :param dict2: dict
    :return: dict
    """

    # format: length_dict = {x: {'ins': 0, 'del': 0} for x in indiv_ids}

    assert sorted(dict1.keys()) == sorted(dict2.keys())

    merge_dict = {x: {'ins': 0, 'del': 0} for x in dict1.keys()}

    for indiv in merge_dict.keys():
        for var in ['ins', 'del']:

            merge_dict[indiv][var] = dict1[indiv][var] + dict2[indiv][var]

    return merge_dict


def sequence_loss_gain(chromo, start, stop, pysam_vcf):

    """
    for a given region returns the number of bases inserted and deleted per individual
    :param chromo: str
    :param start: int
    :param stop: int
    :param pysam_vcf: pysam.VariantFile
    :return: dict
    """

    indiv_ids = pysam_vcf.header.samples
    length_dict = {x: {'ins': 0, 'del': 0} for x in list(indiv_ids) + ['n_indel', 'total_event_length']}

    n_ins, len_ins = 0, 0
    n_del, len_del = 0, 0

    for line in pysam_vcf.fetch(chromo, start, stop):

        # skips individuals with missing genotypes
        if missing_genotypes(line):
            continue

        indel_var = indel_type(line)

        # skip unpolarised sites
        if indel_var is None:
            continue

        allele_lens = {0: len(line.ref) - 1, 1: len(line.alts[0]) - 1}

        # update counters
        if indel_var == 'ins':
            n_ins += 1
            len_ins += abs(len(line.ref) - len(line.alts[0]))
        else:
            n_del += 1
            len_del += abs(len(line.ref) - len(line.alts[0]))

        # get genotype for each indiv
        for indiv in indiv_ids:
            geno = line.samples[indiv]['GT']

            for allele in geno:
                length_dict[indiv][indel_var] += allele_lens[allele]

    # update indel counts
    length_dict['n_indel']['ins'] = n_ins
    length_dict['n_indel']['del'] = n_del
    length_dict['total_event_length']['ins'] = len_ins
    length_dict['total_event_length']['del'] = len_del

    return length_dict


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='vcf file to summarise indel lengths in', required=True)
    parser.add_argument('-call_fa', help='fasta file of callable sites', required=True)
    parser.add_argument('-chr_bed', help='bed file of chromosomes to calc callable sites for', required=True)
    parser.add_argument('-opt_bed', help='Optional bed file of regions to count callables sites, with associated label'
                                         'i.e. /path/to/file.bed.gz,my_sub_region', action='append')
    args = parser.parse_args()

    # file types
    vcf = pysam.VariantFile(args.vcf)
    call_fa = pysam.FastaFile(args.call_fa)
    if args.opt_bed is None:
        region_beds = []
    else:
        region_beds = [(x.split(',')[0], x.split(',')[1]) for x in args.opt_bed]

    # per chromosome
    first = True
    for line in gzip.open(args.chr_bed):
        contig, start_pos, end_pos = line.split()[0], int(line.split()[1]), int(line.split()[2])

        # per chromo genome wide
        seq_lost_gain = sequence_loss_gain(contig, start_pos, end_pos, vcf)
        call_sites = window_call_sites(call_fa, None, [contig, start_pos, end_pos])

        print_dict(seq_lost_gain, contig, 'gwide', call_sites, include_header=first)

        first = False

        # and per chromo per optional regional bed file
        for reg in region_beds:

            first = True
            region = reg[1]
            bed = pysam.TabixFile(reg[0])

            # for each region in bed file
            seq_lost_gain = {}
            call_sites = 0

            for bed_line in bed.fetch(contig, start_pos, end_pos, parser=pysam.asTuple()):

                if first:
                    seq_lost_gain = sequence_loss_gain(bed_line[0], int(bed_line[1]), int(bed_line[2]), vcf)

                else:
                    iterations_loss_gain = sequence_loss_gain(bed_line[0], int(bed_line[1]), int(bed_line[2]), vcf)
                    seq_lost_gain = sum_loss_gain_dict(seq_lost_gain, iterations_loss_gain)

                call_sites += window_call_sites(call_fa, None, [bed_line[0], int(bed_line[1]), int(bed_line[2])])

                first = False

            print_dict(seq_lost_gain, contig, region, call_sites)


if __name__ == '__main__':
    main()
