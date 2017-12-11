#!/usr/bin/env python

import argparse
import os

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', help='Directory containg mafs to be merged', required=True)
parser.add_argument('-out_maf', help='Output directory and maf name', required=True)
args = parser.parse_args()

# variables
maf_source = args.dir
out_maf = args.out_maf
maf_list = sorted([maf_source + maf for maf in os.listdir(maf_source) if maf.endswith('.maf')])

# write output
first_maf = True
with open(out_maf, 'w') as merged_maf:
    for maf in maf_list:
        with open(maf) as maf_in:
            for line in maf_in:
                if line.startswith('#'):
                    if first_maf is True:
                        if line.startswith('##'):
                            merged_maf.write(line)
                else:
                    merged_maf.write(line)
            first_maf = False
            print(maf + ' processed!')

print('Processing complete, merged maf written to ' + out_maf)
