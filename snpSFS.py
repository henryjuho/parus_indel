#!/usr/bin/env python

import argparse
from qsub import *
import sys
import vcf as py_vcf
import re
import numpy
import random

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf',
                    help='VCF file to get site frequency spectrum from',
                    required=True)
parser.add_argument('-folded',
                    help='Specify Y for folded and N for unfolded spectrum',
                    default='Y',
                    choices=['Y', 'N'])
parser.add_argument('-bin',
                    help='Specify genomic regions to bin data by (requires an annotated vcf',
                    default=[],
                    choices=['CDS', 'intron', 'intergenic', 'intergenic_ww_ss', 'zerofold'],
                    action='append')
parser.add_argument('-include_sex', help='If not specified excludes sex chromosomes',
                    default=True, action='store_false')
parser.add_argument('-sfs_out',
                    help='Output sfs data prefix',
                    required=True)
parser.add_argument('-bootstrap', help='Help specify number of times to bootstrap, default is none', default=0)
parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
args = parser.parse_args()

# submission loop
if args.sub is True:
    command_line = [' '.join([x for x in sys.argv if x != '-sub' and x != '-evolgen'])]
    if args.bootstrap is 0:
        time = 8
    else:
        time = 36
    q_sub(command_line, args.sfs_out, t=time, evolgen=args.evolgen)
    sys.exit()


# functions
def bootstrap_sfs(sfs_array, n_boots):
    """
    function that runs a specified number of bootstraps on a sfs and returns the mean and standard error of these

    :param sfs_array:  list [[freq, proportion, bin, normalised_freq, recomb], [...]]
    :param n_boots: integer
    :return : nested list [[freq, mean, bin, normalised_freq, se, norm_se, ci_lwr, ci_upr, ci_lwr_norm, ci_upr_norm],
    [...]]
    """

    region_bin = sfs_array[0][2]

    # create frequency list
    sfs_freqs = []
    for freq in sfs_array:
        indel_total = freq[1]
        for i in range(0, indel_total):
            sfs_freqs.append(freq[0])

    # bootstrap
    bootstrapped_sfs = {f[0]: [] for f in sfs_array}
    for bs in range(0, n_boots):
        # resample with replacement
        resampled_sfs = []
        for i in range(0, len(sfs_freqs)):
            random_no = random.randint(0, len(sfs_freqs)-1)
            resampled_sfs.append(sfs_freqs[random_no])
        for f in bootstrapped_sfs.keys():
            bootstrapped_sfs[f].append(resampled_sfs.count(f))

    # get mean and standard error
    bootstrap_output = []
    # [freq, mean, bin, normalised_freq, se, norm_se, ci_lwr, ci_upr, ci_lwr_norm, ci_upr_norm]
    for f in bootstrapped_sfs.keys():
        mean = numpy.mean(bootstrapped_sfs[f])
        se = numpy.std(bootstrapped_sfs[f])/numpy.sqrt(float(n_boots))
        ci = numpy.percentile(bootstrapped_sfs[f], [2.5, 97.5])

        bootstrap_output.append([f, mean, region_bin, 0, se, 0, ci[0], ci[1], 0, 0])

    return bootstrap_output


# variables
vcf_file = args.vcf
folded = args.folded
output = args.sfs_out
vcf = py_vcf.Reader(open(vcf_file))
number_samples = len(vcf.samples)
bins = args.bin
bootstrap = int(args.bootstrap)
auto_only = args.include_sex
if len(bins) == 0:
    bins.append('no_bins')

# dictionaries
freq_dict = {}

# preload frequency dictionary
# to follow form {bin :{freq: [freq, proportion, bin, normalised_freq]}}
if folded == 'N':
    for region_bin in bins:
        freq_dict[region_bin] = {float(freq)/float(2 * number_samples):
                                 [float(freq)/float(2 * number_samples), 0, region_bin, 0]
                                 for freq in range(1, 2 * number_samples)}
else:
    for region_bin in bins:
        freq_dict[region_bin] = {float(freq)/float(2 * number_samples):
                                 [float(freq)/float(2 * number_samples), 0, region_bin, 0]
                                 for freq in range(1, number_samples + 1)}

# generate folded spectrum if specified
if folded == 'Y':
    for snp in vcf:
        if snp.num_called == number_samples:

            # skip sex chromos
            if auto_only is True:
                if snp.CHROM == 'chrZ' or snp.CHROM == 'chrW':
                    continue

            alt_freq = snp.aaf[0]

            # determine region of variant
            if bins[0] == 'no_bins':
                region = bins[0]
            else:
                try:
                    region = snp.INFO['ANNO']
                except KeyError:
                    continue

            # catch S<->S and W<->W
            alt_seq = str(snp.ALT[0])
            ref_seq = snp.REF
            if 'intergenic_ww_ss' in bins:
                if re.match(r'[CG]', alt_seq) and re.match(r'[CG]', ref_seq) and region == 'intergenic':
                    region = 'intergenic_ww_ss'
                elif re.match(r'[AT]', alt_seq) and re.match(r'[AT]', ref_seq) and region == 'intergenic':
                    region = 'intergenic_ww_ss'

            # catch zerofold
            if region == 'CDS_non_frameshift' and 'zerofold' in bins:
                try:
                    degen = snp.INFO['DEGEN']
                    if degen == 0:
                        region = 'zerofold'
                    else:
                        region = region
                except KeyError:
                    region = region

            # determine minor allele freq
            alt_freq = snp.aaf[0]
            if alt_freq <= 0.5:
                minor_allele_freq = alt_freq
            else:
                minor_allele_freq = 1.0 - alt_freq
            try:
                freq_dict[region][minor_allele_freq][1] += 1
            except KeyError:
                continue

# generate unfolded spectrum if specified
else:
    for snp in vcf:
        if snp.num_called == number_samples:

            # skip sex chromos
            if auto_only is True:
                if snp.CHROM == 'chrZ' or snp.CHROM == 'chrW':
                    continue

            # see if snp is polarised, skip if not
            try:
                ancestral_sequence = snp.INFO['AA']
            except KeyError:
                continue

            alt_seq = str(snp.ALT[0])
            ref_seq = snp.REF
            alt_freq = snp.aaf[0]

            # determine region of variant
            if bins[0] == 'no_bins':
                region = bins[0]
            else:
                try:
                    region = snp.INFO['ANNO']
                except KeyError:
                    continue

            # catch S<->S and W<->W
            if 'intergenic_ww_ss' in bins:
                if re.search(r'[CG]', alt_seq) and re.search(r'[CG]', ref_seq) and region == 'intergenic':
                    region = 'intergenic_ww_ss'
                elif re.search(r'[AT]', alt_seq) and re.search(r'[AT]', ref_seq) and region == 'intergenic':
                    region = 'intergenic_ww_ss'

            # catch zerofold
            if region == 'CDS_non_frameshift' and 'zerofold' in bins:
                try:
                    degen = snp.INFO['DEGEN']
                    if degen == 0:
                        region = 'zerofold'
                    else:
                        region = region
                except KeyError:
                    region = region

            # determine derived allele frequency
            if alt_seq == ancestral_sequence:
                derived_freq = 1.0 - alt_freq
            elif ref_seq == ancestral_sequence:
                derived_freq = alt_freq
            else:
                continue

            # record site frequencies
            try:
                freq_dict[region][derived_freq][1] += 1
            except KeyError:
                continue


# process frequency data - calculate normalised

print '|Folded |Region       | No_SNPs    |\n' \
      '|:------|:-----------:|:----------:|'

new_output = output + '.' + 'folded_' + folded + '_sfs.txt'
with open(new_output, 'w') as sfs:
    if bootstrap is not None:
        sfs.write('\t'.join(['Freq', 'Mean', 'SE', 'Norm', 'Norm_SE', 'Bin',
                             'CI_lwr', 'CI_upr', 'CI_lwr_norm', 'CI_upr_norm']) + '\n')
    else:
        sfs.write('\t'.join(['Freq', 'Proportion', 'Bin', 'Norm']) + '\n')
    for region_key in freq_dict.keys():
        region_dict = freq_dict[region_key]
        total_snps = float(sum([x[1] for x in region_dict.values()]))
        print '| ' + folded + ' |' + region_key + ' |' + str(total_snps) + ' |'
        data = [x for x in region_dict.values()]

        # perform bootstrapping if specified
        if bootstrap != 0:
            bootstrap_data = bootstrap_sfs(data, bootstrap)
            # [freq, mean, bin, normalised_freq, se, norm_se]
            for frequency in bootstrap_data:
                frequency[3] = float(frequency[1])/total_snps
                frequency[5] = frequency[4]/total_snps
                frequency[8] = frequency[6]/total_snps
                frequency[9] = frequency[7]/total_snps
            sorted_data = sorted(bootstrap_data, key=lambda y: y[0])
            for row in sorted_data:
                new_row = [str(z) for z in [row[0], row[1], row[4], row[3], row[5], row[2]] + row[6:]]
                sfs.write('\t'.join(new_row) + '\n')

        # if no bootstrapping specified
        else:
            for frequency in data:
                try:
                    frequency[3] = float(frequency[1])/total_snps
                except ZeroDivisionError:
                    frequency[3] = 0.0
            sorted_data = sorted(data, key=lambda y: y[0])
            for row in sorted_data:
                row = [str(z) for z in row]
                sfs.write('\t'.join(row) + '\n')