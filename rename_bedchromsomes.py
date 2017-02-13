#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-chr_IDs', help='Tab delim text file with two columns in form of chr\trefseqid', required=True)
args = parser.parse_args()

# variables
chr_id_dict = {x.split()[1]: x.split()[0] for x in open(args.chr_IDs)}

# loop through input bed
for line in sys.stdin:
    chromo = line.split()[0]
    try:
        new_chromo = chr_id_dict[chromo]
    except KeyError:
        continue
    new_line = new_chromo + '\t' + '\t'.join(line.split()[1:])
    print(new_line)
