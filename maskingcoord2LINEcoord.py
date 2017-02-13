#!/usr/bin/env python

from __future__ import print_function
import argparse
import sys

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-LINE_type', help='Text file of LINEs to extract coords for', required=True)
args = parser.parse_args()

# variables
lines = {x.rstrip('\n') for x in open(args.LINE_type)}

# loop through stdin
for coord in sys.stdin:
    coord = coord.rstrip('\n').split()
    line_type = coord[3]
    if line_type not in lines:
        continue
    print('\t'.join(coord[0:3]))
