#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import *
import math
import random
import numpy
import sys


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


def read_callable_csv(csv):
    call_sites = {}
    for line in open(csv):
        if not line.startswith('contig'):
            info = line.rstrip().split(',')
            contig, reg, all_call, pol_call = info
            if contig not in call_sites.keys():
                call_sites[contig] = {reg: {'all': float(all_call), 'pol': float(pol_call)}}
            else:
                call_sites[contig][reg] = {'all': float(all_call), 'pol': float(pol_call)}
    return call_sites


def resample_replace(site_freqs):
    resamp_sfs = []
    for i in range(0, len(site_freqs)):
        random_no = random.randint(0, len(site_freqs) - 1)
        resamp_sfs.append(site_freqs[random_no])
    return resamp_sfs


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='Vcf file to get summary stats for', required=True)
    parser.add_argument('-call_csv', help='CSV file of callable sites', required=True)
    parser.add_argument('-mode', help='Variant mode', required=True, choices=['SNP', 'INDEL'])
    parser.add_argument('-bootstrap', help='If specified will perform the requested number of rounds of bootstrapping',
                        default=0)
    parser.add_argument('-no_sex', help='If specified will skip sex chromosomes', default=False, action='store_true')
    parser.add_argument('-per_chromo', help='If specified will print per chromosome stats as well as genome wide',
                        default=False, action='store_true')
    parser.add_argument('-md', help='If specified will output in markdown table format', default=False,
                        action='store_true')
    parser.add_argument('-freq_filter_off', help='Adds multiallelic flag', default=False, action='store_true')
    parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
    parser.add_argument('-evolgen', help='If specified will submit to lab queue', default=False, action='store_true')
    parser.add_argument('-out', help='Output file if submitted')
    args = parser.parse_args()

    # submission loop
    if args.sub is True:
        if args.out is None:
            sys.exit('-out must be specified in conjunction with -sub')
        command_line = [' '.join([x for x in sys.argv if x != '-sub' and x != '-evolgen']) + ' > ' + args.out]
        if args.bootstrap != 0:
            t = 48
        else:
            t = 8
        q_sub(command_line, out=args.out, evolgen=args.evolgen, t=t)
        sys.exit()

    # variables
    vcf_file = args.vcf
    callable_sites = read_callable_csv(args.call_csv)
    mode = args.mode
    markdown = args.md
    bootstrap = int(args.bootstrap)
    no_sex = args.no_sex
    if args.freq_filter_off:
        freq_flag = ' -multi_allelic'
    else:
        freq_flag = ''
    if no_sex:
        sex_flag = '-auto_only'
    else:
        sex_flag = ''
    n = (int(popen_grab('zgrep ^#CHROM {} | wc -w'.format(vcf_file))[0]) - 9) * 2  # 9 columns before sample IDs in VCF

    if args.per_chromo:
        chromo_list = popen_grab('zgrep -v ^# {} | cut -f 1 | uniq'.format(vcf_file))
    else:
        chromo_list = []

    sex_chromos = {'chrZ', 'Z', 'chrW', 'W', 'X', 'XHet', 'Y', 'YHet'}

    if mode == 'SNP':
        regions = ['ALL', 'CDS', 'intron', 'intergenic', 'non-coding', 'AR']  # , 'zerofold', 'fourfold']
    else:
        regions = ['ALL', 'CDS', 'CDS_frameshift', 'CDS_non_frameshift', 'intron', 'intergenic', 'non-coding', 'AR']

    # write header
    if markdown is True:
        print('|region|bin|type|seg_sities|callable|theta_w|t_lwr|t_upr|pi|pi_lwr|pi_upr|tajD|tajD_lwr|tajD_upr|\n'
              '|:----|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|')
    else:
        print('region\tbin\ttype\tseg_sites\tcallable\t'
              'theta_w\tt_lwr\tt_upr\t'
              'pi\tpi_lwr\tpi_upr\t'
              'tajD\ttajD_lwr\ttajD_upr')

    # use sfs script to get freqs per region
    for chromo in chromo_list + ['ALL']:

        # skip sex chromosome if specified
        if no_sex and chromo in sex_chromos:
            continue

        # for each region
        for region in regions:

            if region == 'ALL':
                region_flag = ''
            elif region == 'CDS':
                region_flag = ' -region CDS_frameshift -region CDS_non_frameshift'
            elif region == 'AR':
                region_flag = ' -region intergenic_ar -region intron_ar'
            elif region == 'non-coding':
                region_flag = ' -region intergenic -region intron'
            elif region == 'zerofold':
                region_flag = ' -degen 0'
            elif region == 'fourfold':
                region_flag = ' -degen 4'
            else:
                region_flag = ' -region ' + region

            if region == 'non-coding':
                all_callable = callable_sites[chromo]['intergenic']['all'] + callable_sites[chromo]['intron']['all']
                pol_callable = callable_sites[chromo]['intergenic']['pol'] + callable_sites[chromo]['intron']['pol']
            else:
                all_callable = callable_sites[chromo][region.split('_')[0]]['all']
                pol_callable = callable_sites[chromo][region.split('_')[0]]['pol']

            if mode == 'SNP':
                snp_sfs = popen_grab('~/sfs_utils/vcf2raw_sfs.py -vcf {}{} -chr {}{} -mode snp -folded {}'
                                     .format(vcf_file, freq_flag, chromo, region_flag, sex_flag))

                sfs_list = [([float(x) for x in snp_sfs], 'snp', all_callable)]

            else:
                del_sfs = popen_grab('~/sfs_utils/vcf2raw_sfs.py -vcf {}{} -chr {}{} -mode del {}'
                                     .format(vcf_file, freq_flag, chromo, region_flag, sex_flag))
                ins_sfs = popen_grab('~/sfs_utils/vcf2raw_sfs.py  -vcf {}{} -chr {}{} -mode ins {}'
                                     .format(vcf_file, freq_flag, chromo, region_flag, sex_flag))
                indel_sfs = popen_grab('~/sfs_utils/vcf2raw_sfs.py  -vcf {}{} -chr {}{} -mode indel '
                                       '-folded {}'
                                       .format(vcf_file, freq_flag, chromo, region_flag, sex_flag))

                sfs_list = [([float(x) for x in del_sfs], 'del', pol_callable),
                            ([float(x) for x in ins_sfs], 'ins', pol_callable),
                            ([float(x) for x in indel_sfs], 'indel', all_callable)]

            # process all mute type sfs
            for mute_sfs in sfs_list:

                try:
                    tw = theta_w(n, len(mute_sfs[0])) / mute_sfs[2]
                    pi_val = pi(n, mute_sfs[0]) / mute_sfs[2]
                    tajd = tajimas_d(n, mute_sfs[0])
                except ZeroDivisionError:
                    tw, pi_val, tajd = 0.0, 0.0, 0.0

                # bootstrapping

                if bootstrap == 0 or mute_sfs[2] == 0.0 or tw == 0.0:

                    # set cis for output if no bootstrapping
                    ci_tw = [0, 0]
                    ci_pi = [0, 0]
                    ci_tajd = [0, 0]

                else:
                    bs_theta = []
                    bs_pi = []
                    bs_tajd = []

                    # resample with replacement - calc stats for each resample
                    for bs in range(0, bootstrap):

                        resampled_sfs = resample_replace(mute_sfs[0])

                        bs_theta.append(theta_w(n, len(resampled_sfs)) / mute_sfs[2])
                        bs_pi.append(pi(n, resampled_sfs) / mute_sfs[2])
                        bs_tajd.append(tajimas_d(n, resampled_sfs))

                    # calc cis
                    ci_tw = numpy.percentile(bs_theta, [2.5, 97.5])
                    ci_pi = numpy.percentile(bs_pi, [2.5, 97.5])
                    ci_tajd = numpy.percentile(bs_tajd, [2.5, 97.5])

                if markdown is True:
                    sep = '|'
                    pad = '|'
                    dp = 5
                else:
                    sep = '\t'
                    pad = ''
                    dp = 10
                print(pad + sep.join([chromo, region, mute_sfs[1], str(len(mute_sfs[0])), str(mute_sfs[2]),
                      str(round(tw, dp)), str(round(ci_tw[0], dp)), str(round(ci_tw[1], dp)),
                      str(round(pi_val, dp)), str(round(ci_pi[0], dp)), str(round(ci_pi[1], dp)),
                      str(round(tajd, dp)), str(round(ci_tajd[0], dp)), str(round(ci_tajd[1], dp))]) + pad)

if __name__ == '__main__':
    main()
