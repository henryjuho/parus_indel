#!/usr/bin/env python

import argparse
from hen_utils import *
import os

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-maxL_dir', help='Directory containing maxL.txt files', required=True)
args = parser.parse_args()

# variables
ml_txts = [args.maxL_dir + x for x in os.listdir(args.maxL_dir) if x.endswith('maxL.txt')]

# concat output
counter = 0
for maxL in ml_txts:
    counter += 1
    info = extract_filename(maxL).split('.')
    r = info[-5]
    t = info[-4]
    rep = info[-3]
    with open(maxL) as mltxt:
        for line in mltxt:
            line = line.rstrip('\n')
            if line.startswith('model'):
                if counter == 1:
                    print line + '\tregion\ttest\trep'
            else:
                print line + '\t' + r + '\t' + t + '\t' + rep
