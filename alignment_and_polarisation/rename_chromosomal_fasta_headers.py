#!/usr/bin/env python

import argparse
import os

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fas_in', help='Input fasta directory', required=True)
args = parser.parse_args()

# variables
fastas = [args.fas_in + fa for fa in os.listdir(args.fas_in) if fa.endswith('fa')]

# loop through list
for fasta in fastas:
    new_fasta = open(fasta.replace('.fa', '.rename.fa'), 'w')
    chromo = fasta.split('.')[-2]
    for line in open(fasta):
        if line.startswith('>'):
            new_fasta.write('>' + chromo + '\n')
        else:
            new_fasta.write(line)
    new_fasta.close()
