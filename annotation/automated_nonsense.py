#!/usr/bin/env python

import argparse
import subprocess
from qsub import q_sub


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-cds_fa', help='Fasta file with CDS sequences in', required=True)
    parser.add_argument('-vcf', help='SNP vcf path', required=True)
    parser.add_argument('-out', help='output file stem', required=True)
    parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
    args = parser.parse_args()

    # get chromosome list
    grep_cmd = 'zgrep -v ^# {} | cut -f 1 | uniq'.format(args.vcf)
    chromo_list = subprocess.Popen(grep_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0].split('\n')[:-1]
    chromo_list = [x for x in chromo_list if x.startswith('chr')]

    # loop through chromo list and submit job for each
    for chromo in chromo_list:
        stem = '_'.join([args.out, chromo])

        nonsense_cmd = ('~/parus_indel/annotation/prem_stops_to_bed.py '
                        '-cds_fa {} '
                        '-vcf {} '
                        '-chr {} '
                        '-out {}').format(args.cds_fa, args.vcf, chromo, args.out)

        q_sub([nonsense_cmd],
              out=stem,
              t=48,
              evolgen=args.evolgen)

if __name__ == '__main__':
    main()
