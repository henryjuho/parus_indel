#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import q_sub


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-call_fa', help='Callable fasta file', required=True)
    parser.add_argument('-chr_list', help='File of chromosomes to calc callable sites for', required=True)
    parser.add_argument('-region_list', help='text file with pairs of labels and bed files', required=True)
    parser.add_argument('-out_pre', help='output path and prefix', required=True)
    parser.add_argument('-evolgen', help='If specified will submit to lab queue', default=False, action='store_true')
    args = parser.parse_args()

    for region in open(args.region_list):

        tag, bed = region.split()

        out_stem = args.out_pre + '_' + tag

        cmd = ('~/parus_indel/summary_analyses/callable_sites_bed.py '
               '-call_fa {} -chr_list {} '
               '-bed {} -tag {} '
               '> {}').format(args.call_fa, args.chr_list, bed, tag, out_stem + '.txt')

        q_sub([cmd], out=out_stem, mem=10, rmem=10, evolgen=args.evolgen)


if __name__ == '__main__':
    main()
