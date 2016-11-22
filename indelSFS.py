#!/usr/bin/env python

import argparse
import vcf as py_vcf
import sys
from qsub import *
import random
import numpy

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf',
                    help='VCF file to get site frequency spectrum from',
                    required=True)
parser.add_argument('-folded',
                    help='Specify Y for folded and N for unfolded spectrum',
                    default='Y',
                    choices=['Y', 'N'])
parser.add_argument('-fold_type', help='Type of folded spectra, minior allele or shorter allele',
                    choices=['minor', 'short'], default='short')
parser.add_argument('-bin',
                    help='Specify genomic regions to bin data by (requires an annotated vcf',
                    default=[],
                    choices=['CDS_non_frameshift', 'CDS_frameshift', 'intron', 'intergenic', 'CDS'],
                    action='append')
parser.add_argument('-rbin', help='Type of recombination binning to perform', choices=['None', 'crude', 'poly'],
                    default='None')
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
    if args.bootstrap == 0:
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
    :return : nested list [[freq, mean, bin, normalised_freq, recomb, se, norm_se, ci_lwr, ci_upr, ci_lwr_norm,
    ci_upr_norm], [...]]
    """

    region_bin = sfs_array[0][2]
    recomb_region_bin = sfs_array[0][4]

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
    # [freq, mean, bin, normalised_freq, recomb, se, norm_se, ci_lwr, ci_upr, ci_lwr_norm, ci_upr_norm]
    for f in bootstrapped_sfs.keys():
        mean = numpy.mean(bootstrapped_sfs[f])
        se = numpy.std(bootstrapped_sfs[f])/numpy.sqrt(float(n_boots))
        ci = numpy.percentile(bootstrapped_sfs[f], [2.5, 97.5])

        bootstrap_output.append([f, mean, region_bin, 0, recomb_region_bin, se, 0, ci[0], ci[1], 0, 0])

    return bootstrap_output

# variables
vcf_file = args.vcf
folded = args.folded
folded_type = args.fold_type
output = args.sfs_out
bins = args.bin
auto_only = args.include_sex
bootstrap = int(args.bootstrap)
rr_dict = {}
if len(bins) == 0:
    bins.append('no_bins')
if args.rbin == 'crude':
    rbins = ['micro', 'macro_mid', 'macro_end']
elif args.rbin == 'poly':
    rbins = ['low', 'low_mid', 'mid', 'mid_high', 'high']
    # get variant pos + recomb rate
    indel_rates = [('_'.join([line.CHROM, str(line.POS)]), line.INFO['RR'])
                   for line in py_vcf.Reader(open(vcf_file)) if 'RR' in line.INFO.keys()]
    # sort list
    sorted_rates = sorted(indel_rates, key=lambda w: w[1])
    no_indels_per_bin = len(sorted_rates)/5

    # load rr dict
    n = 0
    for rbin in rbins:
        n_plus_i = n + no_indels_per_bin + 1
        for site in sorted_rates[n:n_plus_i]:
            rr_dict[site[0]] = rbin
        n = n_plus_i

else:
    rbins = ['None']
vcf = py_vcf.Reader(open(vcf_file))
number_samples = len(vcf.samples)

# dictionaries
freq_dict = {}
del_freq_dict = {}
ins_freq_dict = {}

# preload frequency dictionary
# to follow form {recomb: {bin :{freq: [freq, proportion, bin, normalised_freq, recomb]}}}
if folded == 'Y' and folded_type == 'minor':
    no_frequencies = number_samples + 1
else:
    no_frequencies = 2 * number_samples

for recomb in rbins:
    freq_dict[recomb] = {region_bin: {float(freq)/float(2 * number_samples):
                         [float(freq)/float(2 * number_samples), 0, region_bin, 0, recomb]
                         for freq in range(1, no_frequencies)} for region_bin in bins}
    del_freq_dict[recomb] = {region_bin: {float(freq)/float(2 * number_samples):
                             [float(freq)/float(2 * number_samples), 0, region_bin, 0, recomb]
                             for freq in range(1, no_frequencies)} for region_bin in bins}
    ins_freq_dict[recomb] = {region_bin: {float(freq)/float(2 * number_samples):
                             [float(freq)/float(2 * number_samples), 0, region_bin, 0, recomb]
                             for freq in range(1, no_frequencies)} for region_bin in bins}

# generate folded spectrum if specified
if folded == 'Y':
    for indel in vcf:
        if indel.num_called == number_samples:

            # skip sex chromos
            if auto_only is True:
                if indel.CHROM == 'chrZ' or indel.CHROM == 'chrW':
                    continue
            alt_len = len(indel.ALT[0])
            ref_len = len(indel.REF)
            alt_freq = indel.aaf[0]

            # determine recomb
            if rbins[0] == 'None':
                recomb_region = rbins[0]
            elif args.rbin == 'poly':
                try:
                    recomb_region = rr_dict[indel.CHROM + '_' + str(indel.POS)]
                except KeyError:
                    continue
            else:
                try:
                    recomb_region = indel.INFO['RBIN']
                except KeyError:
                    continue

            # determine region of variant
            if bins[0] == 'no_bins':
                region = bins[0]
            else:
                try:
                    region = indel.INFO['ANNO']
                except KeyError:
                    continue

            # determine shorter allele freq
            if folded_type == 'short':
                if alt_len > ref_len:
                    shorter_freq = 1.0 - alt_freq
                else:
                    shorter_freq = alt_freq
                try:
                    freq_dict[recomb_region][region][shorter_freq][1] += 1
                except KeyError:
                    continue

            # determine minor allele freq
            else:
                if alt_freq <= 0.5:
                    minor_allele_freq = alt_freq
                else:
                    minor_allele_freq = 1.0 - alt_freq
                try:
                    if 'CDS' in bins and region.starts('CDS'):
                        freq_dict[recomb_region]['CDS'][minor_allele_freq][1] += 1

                    freq_dict[recomb_region][region][minor_allele_freq][1] += 1
                except KeyError:
                    continue

# generate unfolded spectrum if specified
else:
    for indel in vcf:
        if indel.num_called == number_samples:

            # skip sex chromos
            if auto_only is True:
                if indel.CHROM == 'chrZ' or indel.CHROM == 'chrW':
                    continue

            # see if indel is polarised
            try:
                ancestral_sequence = indel.INFO['AA']
            # skip if unpolarised
            except KeyError:
                continue

            # determine recomb
            if rbins[0] == 'None':
                recomb_region = rbins[0]
            elif args.rbin == 'poly':
                try:
                    recomb_region = rr_dict[indel.CHROM + '_' + str(indel.POS)]
                except KeyError:
                    continue
            else:
                try:
                    recomb_region = indel.INFO['RBIN']
                except KeyError:
                    continue

            alt_seq = indel.ALT[0]
            ref_seq = indel.REF
            alt_freq = indel.aaf[0]
            indel_type = 'Error'

            # determine region of variant
            if bins[0] == 'no_bins':
                region = bins[0]
            else:
                try:
                    region = indel.INFO['ANNO']
                except KeyError:
                    continue

            # determine indel type and derived allele frequency
            if len(alt_seq) == len(ancestral_sequence):
                derived_freq = 1.0 - alt_freq
                if len(alt_seq) > len(ref_seq):
                    indel_type = 'Del'
                elif len(alt_seq) < len(ref_seq):
                    indel_type = 'Ins'
            elif len(ref_seq) == len(ancestral_sequence):
                derived_freq = alt_freq
                if len(alt_seq) < len(ref_seq):
                    indel_type = 'Del'
                elif len(alt_seq) > len(ref_seq):
                    indel_type = 'Ins'

            # catch errors
            if indel_type == 'Error':
                # print indel, ancestral_sequence
                # print alt_seq, len(alt_seq), ref_seq, len(ref_seq), alt_freq
                continue

            # record site frequencies
            if indel_type == 'Del':
                try:
                    if 'CDS' in bins and region.startswith('CDS'):
                        del_freq_dict[recomb_region]['CDS'][derived_freq][1] += 1

                    del_freq_dict[recomb_region][region][derived_freq][1] += 1
                except KeyError:
                    continue
            else:
                try:
                    if 'CDS' in bins and region.startswith('CDS'):
                        ins_freq_dict[recomb_region]['CDS'][derived_freq][1] += 1

                    ins_freq_dict[recomb_region][region][derived_freq][1] += 1
                except KeyError:
                    continue

# process frequency data - calculate normalised
if folded == 'Y':
    dict_list = [['folded', freq_dict]]
else:
    dict_list = [['insertions', ins_freq_dict], ['deletions', del_freq_dict]]

print '|Type   |Recomb  |Region       | No_INDELs  |\n' \
      '|:------|:------:|:-----------:|:----------:|'

for spectrum in dict_list:
    new_output = output + '.' + spectrum[0] + '_sfs.txt'
    with open(new_output, 'w') as sfs:
        if bootstrap != 0:
            sfs.write('\t'.join(['Freq', 'Mean', 'SE', 'Norm', 'Norm_SE', 'Bin', 'Recomb',
                                 'CI_lwr', 'CI_upr', 'CI_lwr_norm', 'CI_upr_norm']) + '\n')
        else:
            sfs.write('\t'.join(['Freq', 'Proportion', 'Bin', 'Norm', 'Recomb']) + '\n')
        for recomb_region_key in spectrum[1].keys():
            recomb_region_dict = spectrum[1][recomb_region_key]
            for region_key in recomb_region_dict.keys():
                region_dict = recomb_region_dict[region_key]
                total_indels = float(sum([x[1] for x in region_dict.values()]))
                print '|' + spectrum[0] + ' |' + recomb_region_key + '|' + region_key + ' |' + str(total_indels) + ' |'
                data = [x for x in region_dict.values()]

                # perform bootstrapping if specified
                if bootstrap != 0:
                    bootstrap_data = bootstrap_sfs(data, bootstrap)
                    # [freq, mean, bin, normalised_freq, recomb, se, norm_se, ci_lwr, ci_upr, ci_lwr_norm, ci_upr_norm]
                    for frequency in bootstrap_data:
                        frequency[3] = float(frequency[1])/total_indels
                        frequency[6] = frequency[5]/total_indels
                        frequency[9] = frequency[7]/total_indels
                        frequency[10] = frequency[8]/total_indels
                    sorted_data = sorted(bootstrap_data, key=lambda y: y[0])
                    for row in sorted_data:
                        new_row = [str(z) for z in [row[0], row[1], row[5], row[3], row[6], row[2], row[4]] + row[7:]]
                        sfs.write('\t'.join(new_row) + '\n')

                # if no bootstrapping specified
                else:
                    for frequency in data:
                        frequency[3] = float(frequency[1])/total_indels
                    sorted_data = sorted(data, key=lambda y: y[0])
                    for row in sorted_data:
                        row = [str(z) for z in row]
                        sfs.write('\t'.join(row) + '\n')
