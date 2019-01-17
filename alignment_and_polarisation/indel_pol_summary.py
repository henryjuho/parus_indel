#!/usr/bin/env python

from __future__ import print_function
import argparse
import subprocess


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='vcf file', required=True)
    parser.add_argument('-align_states', help='alignment states file used in polarisation', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    args = parser.parse_args()

    print('cat', 'count', 'percent', 'region', sep=',')

    # loop through region bed files
    for line in open('regions_beds.txt'):
        region, bed = line.rstrip('\n').split(',')

        # call bedtools intersect to get regional vcf
        vcf_out = '{}{}_indels_tmp.vcf'.format(args.out_dir, region)
        bedtools_cmd = ('bedtools intersect '
                        '-header '
                        '-a {vcf} '
                        '-b {regional_bed} '
                        '> {regional_vcf}').format(vcf=args.vcf, regional_bed=bed, regional_vcf=vcf_out)
        subprocess.call(bedtools_cmd, shell=True)

        # run polarisation script suppressing output vcf
        sum_out = '{}{}_indels_pol_sum.txt'.format(args.out_dir, region)
        pol_cmd = ('~/parus_indel/alignment_and_polarisation/polarise_vcf_indels.py '
                   '-vcf {} '
                   '-align_data {} '
                   '-target_spp Greattit '
                   '-no_vcf > {}').format(vcf_out, args.align_states, sum_out)
        subprocess.call(pol_cmd, shell=True)

        # after each run read in output text files, add region column and percentage column
        total = 0
        for res in open(sum_out):
            if res.startswith('categ'):
                continue

            res = res.split('\t')

            if res[0] == 'total':
                total = float(res[1])

            out_line = [res[0], res[1], (int(res[1])/total) * 100, region]
            print(*out_line, sep=',')


if __name__ == '__main__':
    main()
