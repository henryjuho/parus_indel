#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
import sys
sys.path.append('..')
from summary_analyses.bed_summary_stats import tajimas_d, pi
from vcf2raw_sfs import indel_length, get_derived_freq
from anavar_analyses.sel_vs_neu_anavar import read_callable_txt
from summary_analyses.sfs_correct import correct_sfs


def indel_type(vcf_line):

    """
    polarises indel returns state
    :param vcf_line: pysam.VariantRecord
    :return: str
    """

    ref_seq = vcf_line.ref
    alt_seq = vcf_line.alts[0]
    try:
        anc_seq = vcf_line.info['AA']
    except KeyError:  # ie. not polarised
        return None

    # set derived sequence
    if alt_seq == anc_seq:
        derv_seq = ref_seq
    elif ref_seq == anc_seq:
        derv_seq = alt_seq
    else:
        return None

    # determine type
    if len(anc_seq) > len(derv_seq):  # deletion
        return 'del'
    elif len(anc_seq) < len(derv_seq):  # insertion
        return 'ins'
    else:
        return None


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-call_txt', help=argparse.SUPPRESS,
                        default='../summary_analyses/bgi10_call.txt', required=False)
    parser.add_argument('-region', help='region being analysed, for call sites',
                        choices=['CDS', 'noncoding', 'intergenic', 'introns'], required=True)
    parser.add_argument('-correct_sfs', help='corrects sfs for pol error', default=False, action='store_true')
    args = parser.parse_args()

    # connect to vcf from stdin
    vcf = pysam.VariantFile('-')

    # preload dict
    length_freqs = {x: {i: [] for i in range(1, 51)} for x in ['ins', 'del']}

    # process lines
    for line in vcf:

        # skipZ
        if line.contig == 'chrZ':
            continue

        # get type, derived allele freq and length
        var_type = indel_type(line)

        # skip if not pol
        if var_type is None:
            continue

        length = indel_length(line)
        daf = get_derived_freq(line, var_type, 10)

        # add to dict
        length_freqs[var_type][length].append(daf)

    # get call sites from summary file
    call_sites = read_callable_txt(args.call_txt)['ALL'][args.region]['pol']

    # loop through dict and output length, n, call, pi, tajd
    print('length', 'n_var', 'call', 'pi', 'tajd', 'var_type', sep=',')
    for length in range(1, 51):

        del_freqs = length_freqs['del'][length]
        ins_freqs = length_freqs['ins'][length]

        if args.correct_sfs:
            ins_freqs, del_freqs = correct_sfs(ins_freqs, del_freqs)  # uses default errors as from NC

        for allele_freqs in [[ins_freqs, 'ins'], [del_freqs, 'del']]:

            n_var = len(allele_freqs[0])

            if n_var != 0:
                len_pi = pi(20, allele_freqs[0])/call_sites
                len_d = tajimas_d(20, allele_freqs[0])
            else:
                len_pi = 0
                len_d = 'NA'

            print(length, n_var, call_sites, len_pi, len_d, allele_freqs[1], sep=',')


if __name__ == '__main__':
    main()
