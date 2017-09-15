#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import *
import math


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf',
                        help='Allsites vcf to apply filters to and get callable sites',
                        required=True)
    parser.add_argument('-bed', '--bed_repeats',
                        help='BED file with repeat regions listed',
                        required=True)
    parser.add_argument('-ar_bed', '--ar_bed',
                        help='BED file of ancestral repeats',
                        default='None')
    parser.add_argument('-chr_bed',
                        help='Specifies chromosome to extract callable sites for and coords',
                        default='ALL')
    parser.add_argument('-pol',
                        help='If specified will check if site can be polarised, takes a wga bed file',
                        default='None')
    parser.add_argument('-out_pre',
                        help='output path and prefix',
                        required=True)
    args = parser.parse_args()

    # variables
    all_sites = args.vcf
    repeat_bed = args.bed_repeats
    line_bed = args.ar_bed
    chr_bed = [(x.split()[0], x.split()[1], x.split()[2]) for x in open(args.chr_bed)]
    pol = args.pol
    out_pre = args.out_pre
    file_list = open(out_pre + '_falist.txt', 'w')
    jids = []

    # process each chromosome
    for entry in chr_bed:
        chromo, start, stop = entry[0], int(entry[1]), int(entry[2])

        # min chunk size ~10Mb
        no_chunks = int(math.ceil(stop / 10e6))

        chunks = []
        for i in range(0, no_chunks):
            chunk_size = stop / no_chunks

            if i == no_chunks - 1:
                chunks.append((i, i * chunk_size, stop))

            else:
                chunks.append((i, i * chunk_size, ((i+1) * chunk_size)))

        for part in chunks:

            jid = '{}_part{}'.format(chromo, part[0])
            jids.append(jid)

            fa_out = '{}_{}_part{}.fa'.format(out_pre, chromo, part[0])
            print(fa_out, file=file_list)

            cmd_line = ('~/parus_indel/callable_sites_from_vcf_regional.py '
                        '-vcf {vcf} '
                        '-bed {bed} '
                        '-ar_bed {ar} '
                        '-chr {ch} '
                        '-start {st_pos} '
                        '-end {end_pos} '
                        '-pol {wga} '
                        '> {out}'
                        '').format(vcf=all_sites,
                                   bed=repeat_bed,
                                   ar=line_bed,
                                   ch=chromo,
                                   st_pos=part[1],
                                   end_pos=part[2],
                                   wga=pol,
                                   out=fa_out)
            q_sub([cmd_line], out=out_pre, t=48, evolgen=True, jid=jid)

    file_list.close()

    # concat job
    cat_cmd = 'cat {} | ~/parus_indel/fa_cat.py > {}.fa'.format(file_list, out_pre)
    q_sub([cat_cmd], out=out_pre + '_cat', hold=jids, evolgen=True)

if __name__ == '__main__':
    main()
