#!/usr/bin/env python

import argparse
import gzip

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gff', help='GFF file with annotation data in', required=True)
parser.add_argument('-out', help='Output directory', required=True)
args = parser.parse_args()

# variables
gff_file = args.gff
out_dir = args.out
bed_out = out_dir + 'gene_list.bed'
chromo_ids = {}

# make chromosome dict
for line in gzip.open(gff_file):
    if line.startswith('#'):
        continue
    elif line.split()[2] == 'region':
        try:
            line = line.split()
            ID = line[0]
            chromo = 'chr' + line[8].split(';')[3].strip('chromosome=').strip('Name=')
            chromo_ids[ID] = chromo
        except IndexError:
            continue

# write new gff
with open(bed_out, 'w') as bed:
    for line in gzip.open(gff_file):
        if not line.startswith('#'):
            line = line.split()
            if line[2] == 'gene':
                ID = line[0]
                bed.write(chromo_ids[ID] + '\t' + line[3] + '\t' + line[4] + '\n')

print('Done')
