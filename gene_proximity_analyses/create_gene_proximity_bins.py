#!/usr/bin/env python

from __future__ import print_function
import argparse
import math
import sys
import subprocess


def update_bin_coords(bins, new_bins):

    """
    adds interval bins to overall dict
    :param bins: dict
    :param new_bins: dict
    :return: None
    """

    for bin_id in new_bins.keys():

        if bin_id not in bins.keys():
            bins[bin_id] = new_bins[bin_id]

        else:
            bins[bin_id] += new_bins[bin_id]


def range_to_bins(chromo, start, stop, window_size, max_dist):

    """
    takes a bed interval and splits it into bins of distance
    from each end of the interval, the distance determined
    by the window size
    :param chromo: str
    :param start: int
    :param stop: int
    :param window_size: int
    :param max_dist: int
    :return: dict
    """

    range_len = stop - start

    n_wind = range_len / float(window_size)

    n_bins = math.ceil(n_wind/2.0)

    bin_ranges = {}

    for i in range(0, int(n_bins)):

        bin_id = i + 1

        offset = i * window_size

        # terminate when ben is beyond max dist
        if max_dist is not None:

            if bin_id * window_size > max_dist:
                break

        # when reach midrange and final bin just take remaining middle of range
        if i == n_bins - 1:
            bin_start = start + offset
            bin_end = stop - offset

            bin_ranges[bin_id] = [(chromo, bin_start, bin_end)]

        # for all bins except last (moving towards middle of range from each end)
        else:
            forward_bin_start = start + offset
            forward_bin_end = forward_bin_start + window_size

            backward_bin_end = stop - offset
            backward_bin_start = backward_bin_end - window_size

            bin_ranges[bin_id] = [(chromo, forward_bin_start, forward_bin_end),
                                  (chromo, backward_bin_start, backward_bin_end)]

    return bin_ranges


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-bin_size', help='Size of proximty bins in bp', default=1, type=int)
    parser.add_argument('-out_prefix', help='output location an file prefix', required=True)
    parser.add_argument('-max_distance', help='maximum distance from region desired', type=int)
    args = parser.parse_args()

    bin_size = args.bin_size

    bin_dict = {}

    for line in sys.stdin:

        contig, start, stop = line.split()[0], int(line.split()[1]), int(line.split()[2])

        line_bins = range_to_bins(contig, start, stop, bin_size, args.max_distance)

        update_bin_coords(bin_dict, line_bins)

    # write out bed files and bgzip and index
    for dist in bin_dict.keys():

        dist_out = '{}.bin{}.bed'.format(args.out_prefix, dist)

        with open(dist_out, 'w') as out:

            for coords in bin_dict[dist]:
                print(*coords, sep='\t', file=out)

        subprocess.call('bgzip {}'.format(dist_out), shell=True)
        subprocess.call('tabix -pbed {}'.format(dist_out+'.gz'), shell=True)


if __name__ == '__main__':
    main()
