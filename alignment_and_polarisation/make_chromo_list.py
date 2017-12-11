#!/usr/bin/env python

import argparse
import os
import re

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', help='Input directory containing chromosomal fastas', required=True)
parser.add_argument('-spp', help='Species nickname', required=True)
args = parser.parse_args()

# variables
nickname_prefix = args.spp
input_dir = args.dir

# generate directory contents
chromosomes = sorted([fa for fa in os.listdir(input_dir) if
                     re.search(r'\.chromosome\.[chr]{0,3}([ZLGE]{0,3}[\d]{0,2}[AB]?)\.rename\.fa', fa)])

# write chromosomal list files
with open(input_dir + nickname_prefix + '.chromosome.list', 'w') as output:
    for chromo in chromosomes:
        chromo_no = re.search(r'\.chromosome\.[chr]{0,3}([ZLGE]{0,3}[\d]{0,2}[AB]?)\.rename\.fa', chromo).group(1)
        nickname = nickname_prefix + '.chr' + chromo_no
        output.write(input_dir + chromo + '\t' + nickname + '\n')
