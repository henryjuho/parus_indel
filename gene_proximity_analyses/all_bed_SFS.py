#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='vcf file with indels', required=True)
    parser.add_argument('-bed_list', help='list of bin beds', required=True)
    parser.add_argument('-call_fa', help='callable sites fasta file', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
    args = parser.parse_args()

    cmd = 'zgrep ^{chromo} {bed} | python addSFS.py {vcf} {call_fa} > {output}'

    for dist_bed in open(args.bed_list):

        basename = dist_bed.split('/')[-1].replace('.bed.gz', '')
        dist_bed_pysam = pysam.TabixFile(dist_bed)
        out_stem = '{}{}_sfs'.format(args.out_dir, basename)

        for chromo in dist_bed_pysam.contigs:

            output_name = '{}_{}.bed'.format(out_stem, chromo)

            sfs_cmd = cmd.format(chromo=chromo, bed=dist_bed, vcf=args.vcf, call_fa=args.call_fa, output=output_name)
            print(sfs_cmd)

        # gather cmds
        gather_bed = out_stem + '.bed.gz'
        gather = 'cat {}* | sort -k1,1 -k2,2n | bgzip -c > {}'.format(out_stem, gather_bed)
        tabix = 'tabix -pbed {}'.format(gather_bed)
        print(gather, tabix)


if __name__ == '__main__':
    main()
