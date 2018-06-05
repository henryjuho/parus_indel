#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import q_sub


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-wga', help='Whole genome alignment bed file', required=True)
    parser.add_argument('-region_list', help='Coordinates to calc divergence for, bed format', required=True)
    parser.add_argument('-out_dir', help='Output path and file', required=True)
    parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
    args = parser.parse_args()

    for line in open(args.region_list):

        region, bed = line.split()

        out_stem = '{}gt_indel_div_{}'.format(args.out_dir, region)

        div_cmd = ('~/parus_indel/summary_analyses/indel_divergence.py '
                   '-wga {} -bed {} -tag {} > {}').format(args.wga, bed, region, out_stem + '.txt')

        q_sub([div_cmd], out=out_stem, t=24, evolgen=args.evolgen)

if __name__ == '__main__':
    main()
