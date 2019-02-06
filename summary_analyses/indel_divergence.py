#!/usr/bin/env python

from __future__ import print_function
import argparse
import subprocess
import sys


def popen_grab(cmd):
    output_lines_list = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    return output_lines_list


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-wga', help='Whole genome alignment bed file', required=True)
    parser.add_argument('-bed', help='Coordinates to calc divergence for, bed format', required=True)
    parser.add_argument('-tag', help='region name', required=True)
    args = parser.parse_args()

    # callable sites
    callable_cmd = ('bedtools intersect -a {wga} -b {bed} | '
                    'grep -v ^chrZ | '
                    '~/WGAbed/wga_bed_summary.py -callable'
                    ).format(wga=args.wga, bed=args.bed)

    print(callable_cmd, file=sys.stderr)

    call_sites = [int(x.split('\t')[1]) for x in popen_grab(callable_cmd)]
    n_sites = sum(call_sites)

    # fixed differences
    inframe = list(range(3, 51, 3))
    shift = list(set(range(1, 51)) - set(inframe))

    if args.tag == 'cds_non_frameshift':
        lengths = ','.join([str(x) for x in inframe])
        len_tag = ' -lengths {}'.format(lengths)
    elif args.tag == 'cds_frameshift':
        lengths = ','.join([str(x) for x in shift])
        len_tag = ' -lengths {}'.format(lengths)
    else:
        len_tag = ''

    # indels
    indel_cmd = ('bedtools intersect -a {wga} -b {bed} | '
                 'grep -v ^chrZ | '
                 '~/WGAbed/wga_bed_indels.py -min_coverage 3 -max_length 50 -ref_specific{lt} | '
                 'wc -l').format(wga=args.wga, bed=args.bed, lt=len_tag)

    print(indel_cmd, file=sys.stderr)

    n_indels = int(popen_grab(indel_cmd)[0])

    # ins
    ins_cmd = ('bedtools intersect -a {wga} -b {bed} | '
               'grep -v ^chrZ | '
               '~/WGAbed/wga_bed_indels.py -min_coverage 3 -max_length 50 -ref_specific{lt} | '
               '~/WGAbed/polarise_wga_ref_indels.py -indel_type insertion | '
               'wc -l').format(wga=args.wga, bed=args.bed, lt=len_tag)

    print(ins_cmd, file=sys.stderr)

    n_ins = int(popen_grab(ins_cmd)[0])

    # del
    del_cmd = ('bedtools intersect -a {wga} -b {bed} | '
               'grep -v ^chrZ | '
               '~/WGAbed/wga_bed_indels.py -min_coverage 3 -max_length 50 -ref_specific{lt} | '
               '~/WGAbed/polarise_wga_ref_indels.py -indel_type deletion | '
               'wc -l').format(wga=args.wga, bed=args.bed, lt=len_tag)

    print(del_cmd, file=sys.stderr)
    # print('\n\n')

    n_del = int(popen_grab(del_cmd)[0])

    # process
    print("category", "variation", "seg_sites", "callable", "divergence", sep='\t')

    for var in [['INDEL', n_indels], ['INS', n_ins], ['DEL', n_del]]:
        if n_sites == 0:
            div = 0.0
        else:
            div = float(var[1])/float(n_sites)

        print(args.tag, var[0], var[1], n_sites, div, sep='\t')


if __name__ == '__main__':
    main()
