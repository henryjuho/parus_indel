#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys
import os
sys.path.insert(0, os.getenv('HOME') + '/parus_indel/anavar_analyses')
from gen_gamma_plot_data import calc_sel_bins


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-res_csv', help='results to get binned DFE for', required=True)
    args = parser.parse_args()

    print('cat', 'prop', 'var', 'dist', 'theta', sep='\t')
    for line in open(args.res_csv):

        if line.startswith('run'):
            continue

        else:

            line = line.rstrip().split(',')

            if line[10] == 'neu':
                continue

            theta, scale, shape, var_type, dist_bin = line[3], float(line[4]), float(line[5]), line[8], line[18]

            dfe = calc_sel_bins(shape, scale, nbins=2)

            for cat in dfe:

                print(cat[0], cat[1], var_type, dist_bin, theta, sep='\t')

if __name__ == '__main__':
    main()
