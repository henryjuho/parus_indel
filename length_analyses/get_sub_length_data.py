#!/usr/bin/env python

from __future__ import print_function
import argparse
import subprocess


def wga_indel_lengths(wga_lines):

    counter = 0
    indel_length = 0

    # total seq, no events,
    for line in wga_lines:

        seq_len = len(line.split('\t')[7].split(',')[0]) - 1

        counter += 1
        indel_length += seq_len

    return indel_length, counter


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-windows', help='windows file', required=True)
    parser.add_argument('-wga_bed', help='wga bed file', required=True)
    parser.add_argument('-bed', help='Coordinates to calc stats for, bed format', required=True)
    args = parser.parse_args()

    for window in open(args.windows):

        window = window.rstrip()
        if 'start' in window:
            print(window + '\ttotal_length\tn_indel\tn_callable\tindel_type')
            continue

        window_info = window.rstrip().split()
        chromo, start, stop = window_info[0:3]

        call_sites = 'tabix {} {}:{}-{} | bedtools intersect -a stdin -b {} | wga_bed_summary.py -callable'.format(
            args.wga_bed, chromo, int(start) - 1, stop, args.bed)

        n_sites = int(subprocess.Popen(call_sites, shell=True,
                                       stdout=subprocess.PIPE).communicate()[0].split('\t')[1])

        indel_extract = ('tabix {} {}:{}-{} | bedtools intersect -a stdin -b {} |'
                         'wga_bed_indels.py -max_length 50 -min_coverage 3 -ref_specific | '
                         'polarise_wga_ref_indels.py -indel_type {} ')

        for indel_type in ['insertion', 'deletion']:
            cmd = indel_extract.format(args.wga_bed, chromo, int(start) - 1, stop, args.bed, indel_type)

            variants = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')

            var_len, n_var = wga_indel_lengths(variants)

            print(window + '{}\t{}\t{}\t{}'.format(var_len, n_var, n_sites, indel_type))


if __name__ == '__main__':
    main()
