#!/usr/bin/env python

import argparse
import vcf as py_vcf
import sys
from qsub import *

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
                    choices=['CDS_non_frameshift', 'CDS_frameshift', 'intron', 'intergenic'],
                    action='append')
parser.add_argument('-rbin', help='Type of recombination binning to perform', choices=['None', 'crude'], default='None')
parser.add_argument('-auto_only', help='By default exclude sex chromosomes', default=True, choices=[True, False])
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
auto_only = args.auto_only
if len(bins) == 0:
    bins.append('no_bins')
if args.rbin == 'crude':
    rbins = ['micro', 'macro_mid', 'macro_end']
else:
    rbins = ['None']

# dictionaries
freq_dict = {}
del_freq_dict = {}
ins_freq_dict = {}

# preload frequency dictionary
# to follow form {recomb: {bin :{freq: [freq, proportion, bin, normalised_freq, recomb]}}}
for recomb in rbins:
    freq_dict[recomb] = {region_bin: {float(freq)/float(2 * number_samples):
                         [float(freq)/float(2 * number_samples), 0, region_bin, 0, recomb]
                         for freq in range(1, 2 * number_samples)} for region_bin in bins}
    del_freq_dict[recomb] = {region_bin: {float(freq)/float(2 * number_samples):
                             [float(freq)/float(2 * number_samples), 0, region_bin, 0, recomb]
                             for freq in range(1, 2 * number_samples)} for region_bin in bins}
    ins_freq_dict[recomb] = {region_bin: {float(freq)/float(2 * number_samples):
                             [float(freq)/float(2 * number_samples), 0, region_bin, 0, recomb]
                             for freq in range(1, 2 * number_samples)} for region_bin in bins}

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
            if alt_len > ref_len:
                shorter_freq = 1.0 - alt_freq
            else:
                shorter_freq = alt_freq
            try:
                freq_dict[recomb_region][region][shorter_freq][1] += 1
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
                    del_freq_dict[recomb_region][region][derived_freq][1] += 1
                except KeyError:
                    continue
            else:
                try:
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
        sfs.write('\t'.join(['Freq', 'Proportion', 'Bin', 'Norm', 'Recomb']) + '\n')
        for recomb_region_key in spectrum[1].keys():
            recomb_region_dict = spectrum[1][recomb_region_key]
            for region_key in recomb_region_dict.keys():
                region_dict = recomb_region_dict[region_key]
                total_indels = float(sum([x[1] for x in region_dict.values()]))
                print '|' + spectrum[0] + ' |' + recomb_region_key + '|' + region_key + ' |' + str(total_indels) + ' |'
                data = [x for x in region_dict.values()]
                for frequency in data:
                    frequency[3] = float(frequency[1])/total_indels
                sorted_data = sorted(data, key=lambda y: y[0])
                for row in sorted_data:
                    row = [str(z) for z in row]
                    sfs.write('\t'.join(row) + '\n')
