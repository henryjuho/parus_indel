#!/usr/bin/env python

import argparse
from qsub import *
import sys
import vcf as py_vcf
import re

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
parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
args = parser.parse_args()

# submission loop
if args.sub is True:
    command_line = [' '.join([x for x in sys.argv if x != '-sub' and x != '-evolgen'])]
    q_sub(command_line, args.sfs_out, evolgen=args.evolgen)
    sys.exit()

# variables
vcf_file = args.vcf
folded = args.folded
output = args.sfs_out
vcf = py_vcf.Reader(open(vcf_file))
number_samples = len(vcf.samples)
bins = args.bin
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
    sfs.write('\t'.join(['Freq', 'Proportion', 'Bin', 'Norm']) + '\n')
    for region_key in freq_dict.keys():
        region_dict = freq_dict[region_key]
        total_snps = float(sum([x[1] for x in region_dict.values()]))
        print '| ' + folded + ' |' + region_key + ' |' + str(total_snps) + ' |'
        data = [x for x in region_dict.values()]
        for frequency in data:
            try:
                frequency[3] = float(frequency[1])/total_snps
            except ZeroDivisionError:
                frequency[3] = 0.0
        sorted_data = sorted(data, key=lambda y: y[0])
        for row in sorted_data:
            row = [str(z) for z in row]
            sfs.write('\t'.join(row) + '\n')
