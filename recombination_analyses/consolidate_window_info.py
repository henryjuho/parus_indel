#!/usr/bin/env python

from __future__ import print_function
from vcf2raw_sfs import vcf2sfs
import gzip
from window_sel_vs_neu_anavar import window_call_sites
import pysam
import sys
sys.path.insert(0, '/Users/henryjuho/parus_indel/summary_analyses')
from summary_stats_gt import theta_w, pi, tajimas_d


def rec_rates(rec_file):

    rec_dict = {}

    for line in open(rec_file):
        if not line.startswith('window'):
            line = line.split()

            rec_dict[line[0]] = line[1]

    return rec_dict


def main():

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
        ins = list(ins)

        dels = vcf2sfs(vcf, mode='del',
                       chromo=chromo, start=int(start), stop=int(stop),
                       regions=['intron', 'intergenic'])
        dels = list(dels)

        indels = vcf2sfs(vcf, mode='indel',
                         chromo=chromo, start=int(start), stop=int(stop),
                         regions=['intron', 'intergenic'])

        n_indels = len(list(indels))
        pol_success = (len(ins) + len(dels)) / float(n_indels)

        # callsites
        n_call = window_call_sites(call_sites, nc_bed, (chromo, int(start), int(stop)))

        if len(ins) == 0 or len(dels) == 0 or n_call == 0:
            ins_t, ins_pi, ins_taj, dels_t, dels_pi, dels_taj = 'NA', 'NA', 'NA', 'NA', 'NA', 'NA'
        else:
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
