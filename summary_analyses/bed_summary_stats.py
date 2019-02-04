#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import *
import math
import sys
sys.path.extend('..')
from sfs_correct import correct_sfs


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
    parser.add_argument('-bed', help='regions to summarise', required=True)
    parser.add_argument('-tag', help='region name', required=True)
    parser.add_argument('-correct_sfs', help='Corrects sfs for pol error', default=False, action='store_true')
    parser.add_argument('-no_z', help=argparse.SUPPRESS, required=False, default=True, action='store_false')
    args = parser.parse_args()

    # variables
    no_sex = args.no_z

    if no_sex:
        sex_flag = ' -auto_only'
    else:
        sex_flag = ''
    n = (int(popen_grab('zgrep ^#CHROM {} | wc -w'.format(args.indel_vcf))[0]) - 9) * 2  # 9 columns before sample IDs

    # pol error dict
    cds_e = (0.0799355752829, 0.0367725125655)
    nc_e = (0.0110086484429, 0.0166354937984)
    ar_e = (0.0302288845982, 0.0261469211997)

    e_dict = {'ALL': nc_e, 'noncoding': nc_e, 'noncoding_noUCNEs': nc_e, 'intergenic': nc_e,
              'introns': nc_e, 'CDS': cds_e, 'cds_frameshift': cds_e, 'cds_non_frameshift': cds_e,
              '0fold': cds_e, '4fold': cds_e, 'nonsense': cds_e, 'UCNE': nc_e, 'AR': ar_e}

    # write header
    print('category', 'variation', 'seg_sites', 'theta_w', 'pi', 'tajD', sep='\t')

    modes = ['snp', 'ins', 'del', 'indel']
    sfs_dict = {x: [] for x in modes}
    for mode in modes:

        if mode == 'snp' or mode == 'indel':
            folded = ' -folded'
        else:
            folded = ''

        if mode == 'snp':
            vcf_file = args.snp_vcf
        else:
            vcf_file = args.indel_vcf

        # insure none slip through in wrong cds
        if args.tag == 'cds_frameshift':
            region = ' -region ' + args.tag
        elif args.tag == 'cds_non_frameshift':
            region = ' -region ' + args.tag
        else:
            region = ''

        sfs_cmd = ('bedtools intersect -header -u -a {vcf} -b {bed} | '
                   '~/sfs_utils/vcf2raw_sfs.py -mode {mode}{fold_flag}{sex_flag}{region_flag}'
                   .format(vcf=vcf_file, bed=args.bed, mode=mode, fold_flag=folded,
                           sex_flag=sex_flag, region_flag=region))

        sfs = [float(x) for x in popen_grab(sfs_cmd)]

        sfs_dict[mode] = sfs

        # correct sfs
        if args.correct_sfs:

            ei, ed = e_dict[args.tag]

            new_ins, new_del = correct_sfs(sfs_dict['ins'], sfs_dict['del'], e_i=ei, e_d=ed)

            sfs_dict['ins'] = new_ins
            sfs_dict['del'] = new_del

    # calc stats
    for var_type in modes:

        var_sfs = sfs_dict[var_type]

        try:
            tw = theta_w(n, len(var_sfs))
            pi_val = pi(n, var_sfs)
            tajd = tajimas_d(n, var_sfs)
        except ZeroDivisionError:
            tw, pi_val, tajd = 0, 0, float('nan')

        print(args.tag, var_type.upper(), len(var_sfs), tw, pi_val, tajd, sep='\t')


if __name__ == '__main__':
    main()
