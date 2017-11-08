#!/usr/bin/env python

from __future__ import print_function
import argparse


def window_dict(wind_file):

    """
    transforms window data to dict
    :param wind_file: file object
    :return: list(str, dict)
    """

    wind_data = wind_file.readlines()
    header = wind_data[0].rstrip()
    wind_dict = {x.split()[3]: x.split() for x in wind_data[1:]}

    return wind_dict, header


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-results', help='anavar results csv file', required=True)
    parser.add_argument('-windows', help='text file of window data', required=True)
    args = parser.parse_args()

    windows = window_dict(open(args.windows))

    for line in open(args.results):
        line = line.rstrip()

        if line.startswith('run'):
            print(*line.split(',') + [windows[1].split('window\t')[1]], sep=',')
        else:
            line = line.split(',')
            line_window = line[-1]
            window_stats = windows[0][line_window][4:]

            print(*line + window_stats, sep=',')


if __name__ == '__main__':
    main()
