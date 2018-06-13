#!/usr/bin/env python

from __future__ import print_function
import argparse
import subprocess


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-windows', help='windows file', required=True)
    parser.add_argument('-wga_bed', help='wga bed file', required=True)
    args = parser.parse_args()

    for window in open(args.windows):

        window = window.rstrip()

        if 'start' in window:
            print(window + '\tn_ins_sub\tn_del_sub')
            continue

        window_info = window.rstrip().split()
        chromo, start, stop = window_info[0:3]

        indel_extract = ('tabix {} {}:{}-{} | '
                         'wga_bed_indels.py -max_length 50 -min_coverage 3 -ref_specific | '
                         'polarise_wga_ref_indels.py -indel_type {} | '
                         'wc -l')

        for indel_type in ['insertion', 'deletion']:
            cmd = indel_extract.format(args.wga_bed, chromo, int(start)-1, stop, indel_type)
            no_of_var = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[0]

            window += '\t' + no_of_var

        print(window)


if __name__ == '__main__':
    main()
