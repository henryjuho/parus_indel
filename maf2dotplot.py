#!/usr/bin/env python

import argparse
import os

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-maf', help='maf file to plot alignments for', required=True)
parser.add_argument('-out', help='output directory', required=True)
args = parser.parse_args()

# variables
maf_in = args.maf
out_dir = args.out
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
out = open(out_dir + maf_in[maf_in.rfind('/')+1:].rstrip('.maf') + '.plot_data', 'w')

# read maf
spp_1 = True
out_line = ''
with open(maf_in) as maf:
    for line in maf:
        if line.startswith('s'):
            line = line.split()
            spp = line[1]
            pos = line[2]
            if spp_1 is True:
                out_line = spp + '\t' + pos + '\t'
            if spp_1 is False:
                out_line += spp + '\t' + pos + '\n'
                out.write(out_line)
            spp_1 = not spp_1

out.close()
