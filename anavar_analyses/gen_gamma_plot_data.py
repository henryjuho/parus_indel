#!/usr/bin/env python

from __future__ import print_function
import argparse
from numpy.random import gamma


def calc_sel_bins(sh, sc, nbins=4):

    total_sim = 100000
    i_gamma = list(gamma(sh, sc, total_sim))
    if nbins == 4:
        sbins = [('0 - 1', -0.0000000000001, 1), ('1 - 10', 1, 10), ('10 - 100', 10, 100), ('>100', 100, 100000000000)]
    else:
        sbins = [('0 - 1', -0.0000000000001, 1), ('>1', 1, 100000000000)]
    # output proportion in each cat
    for sbin in sbins:
        subset = [x for x in i_gamma if sbin[1] < x <= sbin[2]]

        n_in_cat = len(subset)
        proportion = n_in_cat / float(total_sim)

        yield sbin[0], proportion


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-sh_d', help='del shape', type=float, required=True)
    parser.add_argument('-sc_d', help='del scale', type=float, required=True)
    parser.add_argument('-sh_i', help='ins shape', type=float, required=True)
    parser.add_argument('-sc_i', help='ins scale', type=float, required=True)
    args = parser.parse_args()

    print('cat', 'prop', 'var', sep='\t')

    for y in [('DEL', args.sh_d, args.sc_d), ('INS', args.sh_i, args.sc_i)]:
        gamma_dist = calc_sel_bins(y[1], y[2])

        for x in gamma_dist:
            print(x[0], x[1], y[0], sep='\t')

if __name__ == '__main__':
    main()
