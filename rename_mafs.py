#!/usr/bin/env python

import argparse
import os
from shutil import copyfile

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', help='Directory of mafs to rename', required=True)
args = parser.parse_args()

# variables
in_dir = args.dir
mafs = [in_dir + maf for maf in os.listdir(in_dir) if maf.endswith('maf')]
out_dir = in_dir + 'renamed_mafs/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

# rename mafs
for maf in mafs:
    new_maf = out_dir + maf.replace('_vs_', '.')
    copyfile(maf, new_maf)
