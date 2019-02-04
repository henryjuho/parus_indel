#!/usr/bin/env python

from __future__ import print_function
from vcf2raw_sfs import vcf2sfs
import argparse
import gzip
import pysam
import sys
sys.path.append('..')
from summary_analyses.summary_stats_gt import theta_w, pi, tajimas_d
from summary_analyses.sfs_correct import correct_sfs


def window_call_sites(call_fa, region_bed, window_coords):

    """
    returns number of callable sites for specied region in window
    :param call_fa: pysam.FastaFile()
    :param region_bed: pysam.TabixFile()
    :param window_coords: tuple
    :return: int
    """

    if region_bed is None:
        regions = [(window_coords[0], window_coords[1], window_coords[2])]
    else:
        regions = region_bed.fetch(window_coords[0], window_coords[1], window_coords[2], parser=pysam.asTuple())

    call_sites = 0

    for reg in regions:
        call_seq = call_fa.fetch(reg[0], int(reg[1]), int(reg[2]))
        call_sites += call_seq.count('K')

    return call_sites


def rec_rates(rec_file):

    rec_dict = {}

    for line in open(rec_file):
        if not line.startswith('window'):
            line = line.split()

            rec_dict[line[0]] = line[1]

    return rec_dict


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-correct_sfs', default=False, action='store_true', help=argparse.SUPPRESS)
    args = parser.parse_args()

    # paths
    vcf = '/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Analysis_ready_data/' \
          'final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz'

    rec = '/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/anavar_analysis/' \
          'window_analysis/gt_2Mb_window_rec_rates.txt'

    windows = '/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/anavar_analysis/' \
              'window_analysis/gt_windows.2Mb.bed.gz'

    call_sites = pysam.FastaFile('/Users/henryjuho/sharc_fastdata/GT_ref/bgi_10birds.callable.fa')

    noncoding_bed = '/Users/henryjuho/sharc_fastdata/GT_ref/gt_noncoding.bed.gz'
    rec_data = rec_rates(rec)
    nc_bed = pysam.TabixFile(noncoding_bed)

    # header
    print('chr', 'start', 'stop', 'window', 'rec_rate', 'callable',
          'n_ins', 'theta_ins', 'pi_ins', 'tajd_ins',
          'n_del', 'theta_del', 'pi_del', 'tajd_del', 'n_indel', 'pol_success', sep='\t')

    for line in gzip.open(windows):

        chromo, start, stop, wind_id = line.split()

        # sfs
        ins = vcf2sfs(vcf, mode='ins',
                      chromo=chromo, start=int(start), stop=int(stop),
                      regions=['intron', 'intergenic'])

        dels = vcf2sfs(vcf, mode='del',
                       chromo=chromo, start=int(start), stop=int(stop),
                       regions=['intron', 'intergenic'])

        # correct if specified
        if args.correct_sfs:
            ins, dels = correct_sfs(list(ins), list(dels))
        else:
            ins = list(ins)
            dels = list(dels)

        indels = vcf2sfs(vcf, mode='indel', fold=True,
                         chromo=chromo, start=int(start), stop=int(stop),
                         regions=['intron', 'intergenic'])

        n_indels = len(list(indels))

        # callsites
        n_call = window_call_sites(call_sites, nc_bed, (chromo, int(start), int(stop)))

        if len(ins) == 0 or len(dels) == 0 or n_call == 0:
            ins_t, ins_pi, ins_taj, dels_t, dels_pi, dels_taj, pol_success = 'NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
        else:
            pol_success = (len(ins) + len(dels)) / float(n_indels)

            # summary stats
            ins_t = theta_w(20, len(ins)) / float(n_call)
            ins_pi = pi(20, ins) / float(n_call)
            ins_taj = tajimas_d(20, ins)

            dels_t = theta_w(20, len(dels)) / float(n_call)
            dels_pi = pi(20, dels) / float(n_call)
            dels_taj = tajimas_d(20, dels)

        # rec rate
        rr = rec_data[wind_id]

        print(chromo, start, stop, wind_id, rr, n_call,
              len(ins), ins_t, ins_pi, ins_taj,
              len(dels), dels_t, dels_pi, dels_taj,
              n_indels, pol_success, sep='\t')


if __name__ == '__main__':
    main()
