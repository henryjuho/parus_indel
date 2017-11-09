#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
from vcf2raw_sfs import is_auto


def indel_length(vcf_line):

    """
    takes pysam variant and returns length
    :param vcf_line: pysam variant
    :return: int
    """
    ref_len = len(vcf_line.ref)
    alt_len = len(vcf_line.alts[0])

    length = abs(ref_len - alt_len)

    return length


def indel_type(vcf_line):

    """
    function takes a pysam vcf variant and returns ins of del
    :param vcf_line: pysam variant
    :return: str
    """

    ref_seq = vcf_line.ref
    alt_seq = vcf_line.alts[0]

    try:
        anc_seq = vcf_line.info['AA']
    except KeyError:  # ie. not polarised
        return None

    # set derived sequence and freq
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

    else:  # snp
        return None


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='indel only vcf file', required=True)
    parser.add_argument('-region', help='regional results to include', choices=['CDS', 'non-coding'], action='append')
    parser.add_argument('-auto_only', help='If specified restricts to autosomes', default=False, action='store_true')
    args = parser.parse_args()

    # variables
    vcf = pysam.VariantFile(args.vcf)
    regions = ['gwide']
    if args.region is not None:
        regions += args.region

    length_data = {}

    for x in regions:
        length_data[x] = {'indel': {}, 'ins': {}, 'del': {}}

    for indel in vcf.fetch():

        # skips sex chromosomes if specified
        if args.auto_only and not is_auto(indel):
            continue

        ins_or_del = indel_type(indel)
        indel_len = indel_length(indel)

        try:
            indel_region = indel.info['ANNO']
        except KeyError:
            indel_region = None

        for reg in regions:

            # skips cds dict if indel not in cds
            if reg == 'CDS':
                if indel_region != 'CDS_frameshift' and indel_region != 'CDS_non_frameshift':
                    continue

            # skips non-coding if indel not non-coding
            if reg == 'non-coding':
                if indel_region != 'intergenic' and indel_region != 'intron':
                    continue

            reg_dict = length_data[reg]

            all_data = reg_dict['indel']

            # update count of len x indels
            if indel_len not in all_data.keys():
                all_data[indel_len] = 0

            all_data[indel_len] += 1

            # if polarisable update relevant length count
            if ins_or_del is not None:
                pol_dict = reg_dict[ins_or_del]

                if indel_len not in pol_dict.keys():
                    pol_dict[indel_len] = 0

                pol_dict[indel_len] += 1

    # output table
    print('length', 'count', 'variant', 'region', sep='\t')
    for r in regions:
        for t in ['indel', 'ins', 'del']:
            data = length_data[r][t]
            for length in sorted(data.keys()):
                count = data[length]
                print(length, count, t, r, sep='\t')


if __name__ == '__main__':
    main()
