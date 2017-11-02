#!/usr/bin/env python

from __future__ import print_function
from vcf2raw_sfs import vcf2sfs
import gzip


def rec_rates(rec_file):

    rec_dict = {}

    for line in open(rec_file):
        if not line.startswith('window'):
            line = line.split()

            rec_dict[line[0]] = line[1]

    return rec_dict


def main():

    vcf = '/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Analysis_ready_data/' \
          'final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz'

    rec = '/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/anavar_analysis/' \
          'window_analysis/gt_2Mb_window_rec_rates.txt'

    windows = '/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/anavar_analysis/' \
              'window_analysis/gt_windows.2Mb.bed.gz'

    rec_data = rec_rates(rec)

    print('chr', 'start', 'stop', 'window', 'rec_rate', 'n_indels', 'n_ins', 'n_del', sep='\t')

    for line in gzip.open(windows):

        chromo, start, stop, id = line.split()

        # no variants
        indels = vcf2sfs(vcf, mode='indel', fold=True,
                         chromo=chromo, start=int(start), stop=int(stop),
                         regions=['intron', 'intergenic'])
        no_indels = len(list(indels))

        ins = vcf2sfs(vcf, mode='ins',
                      chromo=chromo, start=int(start), stop=int(stop),
                      regions=['intron', 'intergenic'])
        no_ins = len(list(ins))

        dels = vcf2sfs(vcf, mode='del',
                       chromo=chromo, start=int(start), stop=int(stop),
                       regions=['intron', 'intergenic'])
        no_dels = len(list(dels))

        # rec rate
        rr = rec_data[id]

        print(chromo, start, stop, id, rr, no_indels, no_ins, no_dels, sep='\t')


if __name__ == '__main__':
    main()
