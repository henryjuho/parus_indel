#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import q_sub


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-indel_vcf', help='Vcf file to get summary stats for', required=True)
    parser.add_argument('-snp_vcf', help='Vcf file to get summary stats for', required=True)
    parser.add_argument('-region_list', help='text file with pairs of labels and bed files', required=True)
    parser.add_argument('-out_pre', help='output path and prefix', required=True)
    parser.add_argument('-evolgen', help='If specified will submit to lab queue', default=False, action='store_true')
    args = parser.parse_args()

    for region in open(args.region_list):

        tag, bed = region.split()

        out_stem = args.out_pre + '_' + tag

        cmd = ('~/parus_indel/summary_analyses/bed_summary_stats.py '
               '-indel_vcf {} -snp_vcf {} '
               '-bed {} -tag {} '
               '> {}').format(args.indel_vcf, args.snp_vcf, bed, tag, out_stem + '_stats.txt')

        q_sub([cmd], out=out_stem)


if __name__ == '__main__':
    main()
