#!/usr/bin/env python

import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-maf', help='Alignment file to produce summary stats for', required=True)
parser.add_argument('-fa_list', help='List file in format [fasta] tab [nickname_in_maf]', required=True)
args = parser.parse_args()

# variables
maf = args.maf
fastas = [(fa.split()[0], fa.split()[1]) for fa in open(args.fa_list)]

# Calc % aligned
# 1) Get genome lengths
lengths = {}

for fa in fastas:
    nick = fa[1]
    fa = fa[0]
    faCount_out = subprocess.Popen('faCount < ' + fa, shell=True, stdout=subprocess.PIPE)
    sequence_info = (faCount_out.communicate()[0]).split('\n')
    genome_len = [length.split()[1] for length in sequence_info if length.startswith('total')]
    lengths[nick] = genome_len[0]

# 2) Count bp for each spp in alignment
mapped_lengths = {}
with open(maf) as maf_data:
    for line in maf_data:
        if line.startswith('s'):
            line = line.split()
            ID = line[1].split('.')[0].split('_')[0]
            seq_len = len(line[6])
            if ID not in mapped_lengths.keys():
                mapped_lengths[ID] = seq_len
            else:
                mapped_lengths[ID] += seq_len

# 3) Print % aligned
output = ('MappingSummaryData\n\n'
          '------------------------------------------------------\n'
          'Species\tGenome_len\tMapped_len\tPercent_mapped\n'
          '------------------------------------------------------\n')

for spp in lengths.keys():
    gen_bp = lengths[spp]
    mapped_bp = mapped_lengths[spp]
    percent_mapped = int((float(mapped_bp) / float(gen_bp)) * 100.0)
    output += spp + '\t' + str(gen_bp) + '\t' + str(mapped_bp) + '\t' + str(percent_mapped) + '\n'

output += '------------------------------------------------------'

print(output)
