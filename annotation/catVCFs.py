#!/usr/bin/env python

import argparse
import os

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf',
                    help='VCF file to add to list to concatenate, specify multiple, and in order',
                    required=True,
                    action='append')
parser.add_argument('-out_vcf', help='Path and name of output vcf', required=True)
parser.add_argument('-clean', help='If specified will remove unmerged vcfs', action='store_true', default=False)
args = parser.parse_args()

# variables
vcfs = args.vcf
out = open(args.out_vcf, 'w')
clean = args.clean

# loop through vcf
vcf_counter = 0
for vcf in vcfs:
    vcf_counter += 1
    for line in open(vcf):
        if line.startswith('#'):
            if vcf_counter == 1:
                out.write(line)
            else:
                continue
        else:
            out.write(line)

out.close()

# clean up unmerged vcfs
if clean is True:
    for vcf in vcfs:
        os.remove(vcf)
