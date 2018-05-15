#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import *
import math


def popen_grab(cmd):
    output_lines_list = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    return output_lines_list


def theta_w(n, segsites):
    # theta = S/a
    alpha = sum(1.0/z for z in range(1, n))
    theta = float(segsites) / alpha
    return theta


def pi(n, allele_frq_list):
    no_seqs_dif = [(1.0 - raf**2 - (1.0-raf)**2) * (n/(n-1.0)) for raf in allele_frq_list]
    seq_pi = sum(no_seqs_dif)
    return seq_pi


def tajimas_d(n, allele_frq_list):
    segsites = float(len(allele_frq_list))
    little_d = pi(n, allele_frq_list) - theta_w(n, segsites)

    a1 = sum(1.0/z for z in range(1, n))
    a2 = sum(1.0/z**2 for z in range(1, n))

    e1 = (1.0 / a1) * (((n + 1.0) / (3.0 * (n - 1.0))) - (1.0 / a1))
    e2 = (1.0 / (a1**2 + a2)) * \
         (((2.0 * (n**2 + n + 3.0)) / ((9.0 * n) * (n - 1.0))) -
          ((n + 2.0) / (n * a1)) +
          (a2 / a1**2))
    # print(e2)
    vd = (e1 * segsites) + ((e2 * segsites) * (segsites - 1.0))
    big_d = little_d / math.sqrt(vd)

    return big_d


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-indel_vcf', help='Vcf file to get summary stats for', required=True)
    parser.add_argument('-snp_vcf', help='Vcf file to get summary stats for', required=True)
    # parser.add_argument('-region_list', help='text file with pairs of labels and bed files', required=True)
    parser.add_argument('-bed', help='regions to summarise', required=True)
    parser.add_argument('-tag', help='region name', required=True)
    parser.add_argument('-no_z', help=argparse.SUPPRESS, required=False, default=True, action='store_false')
    args = parser.parse_args()

    # variables
    no_sex = args.no_z

    if no_sex:
        sex_flag = ' -auto_only'
    else:
        sex_flag = ''
    n = (int(popen_grab('zgrep ^#CHROM {} | wc -w'.format(args.indel_vcf))[0]) - 9) * 2  # 9 columns before sample IDs

    # write header
    print('category', 'variation', 'seg_sites', 'theta_w', 'pi', 'tajD', sep='\t')

    for mode in ['snp', 'ins', 'del', 'indel']:

        if mode == 'snp' or mode == 'indel':
            folded = ' -folded'
        else:
            folded = ''

        if mode == 'snp':
            vcf_file = args.snp_vcf
        else:
            vcf_file = args.indel_vcf

        sfs_cmd = ('bedtools intersect -a {vcf} -b {bed} | '
                   '~/sfs_utils/vcf2raw_sfs.py -mode {mode}{fold_flag}{sex_flag}'
                   .format(vcf=vcf_file, bed=args.bed, mode=mode, fold_flag=folded, sex_flag=sex_flag))

        print(sfs_cmd)

        sfs = [float(x) for x in popen_grab(sfs_cmd)]

        tw = theta_w(n, len(sfs))
        pi_val = pi(n, sfs)
        tajd = tajimas_d(n, sfs)

        print(args.tag, mode.upper(), len(sfs), tw, pi_val, tajd, sep='\t')


if __name__ == '__main__':
    main()
