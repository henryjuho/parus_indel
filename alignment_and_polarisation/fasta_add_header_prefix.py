#!/usr/bin/env python

import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fa', help='Fasta file to add prefix to headers', required=True)
parser.add_argument('-pre', help='Prefix', required=True)
args = parser.parse_args()

# variables
fa = args.fa
prefix = args.pre
out = fa.rstrip('.fa') + '.rename.fa'

# rename
with open(out, 'w') as new_fa:
    for line in open(fa):
        if line.startswith('>'):
            line = line.replace('>', '>' + prefix)
            new_fa.write(line)
        else:
            new_fa.write(line)

print('Done')
