#!/usr/bin/env python

from __future__ import print_function
import argparse
import subprocess
from qsub import q_sub
import sys


def popen_grab(cmd):
    output_lines_list = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    return output_lines_list


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-wga', help='Whole genome alignment bed file', required=True)
    parser.add_argument('-bed', help='Coordinates to calc divergence for, bed format', required=True)
    parser.add_argument('-chromo', help=argparse.SUPPRESS, default='all')
    parser.add_argument('-out', help='Output path and file', required=True)
    args = parser.parse_args()

    if args.chromo == 'all':
        chromo_list = popen_grab('zcat {} | cut -f 1 | uniq'.format(args.wga))

        with open(args.out, 'w') as out_file:
            print('chromo\tindels\tcallable\tdivergence', file=out_file)

        for x in chromo_list:

            cmd = ' '.join(sys.argv) + ' -chromo {}'.format(x)
            outs = args.out.replace('.txt', '') + '_{}'.format(x)

            q_sub([cmd], out=outs, evolgen=True)

        sys.exit()

    else:
        callable_cmd = 'zgrep -w ^{} {} | bedtools intersect -a stdin -b {} | wga_bed_summary.py -callable'.format(
            args.chromo, args.wga, args.bed)

        print(callable_cmd, file=sys.stdout)

        call_sites = popen_grab(callable_cmd)[0].split('\t')

        n_sites = int(call_sites[1])

        indel_cmd = ('zgrep -w ^{} {} | bedtools intersect -a stdin -b {} | '
                     'wga_bed_indels.py -min_coverage 3 -max_length 50 -ref_specific | wc -l').format(
            args.chromo, args.wga, args.bed)

        print(indel_cmd, file=sys.stdout)

        n_indels = int(popen_grab(indel_cmd)[0])

        if n_sites == 0:
            div = 0.0
        else:
            div = float(n_indels)/float(n_sites)

        with open(args.out, 'a') as out_file:
            print(args.chromo, n_indels, n_sites, div, file=out_file, sep='\t')


if __name__ == '__main__':
    main()
