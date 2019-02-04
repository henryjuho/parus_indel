#!/usr/bin/env python

from __future__ import print_function
import pandas as pd
import sys


def main():

    windows = sys.argv[1]

    table = pd.read_table(windows)
    filtered_table = table[table.n_ins + table.n_del >= 500]
    filtered_table = filtered_table[filtered_table.chr != 'chrZ']

    filtered_table.to_csv(path_or_buf='filtered_' + windows, sep='\t', index=False)


if __name__ == '__main__':
    main()
