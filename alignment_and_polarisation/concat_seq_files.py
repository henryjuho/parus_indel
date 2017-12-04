#!/usr/bin/env python

import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-seq', help='Seq files to concatonate, specify multiple', action='append', required=True)
parser.add_argument('-out', help='Full output path and file name', required=True)
args = parser.parse_args()

# variables
seqs = args.seq
out_name = args.out

# concatonated
with open(out_name, 'w') as out:
    for seq_file in seqs:
        for line in open(seq_file):
            if line.startswith('ID'):
                if seq_file == seqs[0]:
                    out.write(line)
                else:
                    continue
            else:
                out.write(line)
