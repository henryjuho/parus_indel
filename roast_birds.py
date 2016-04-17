#!/usr/bin/env python

import argparse
import os
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-maf_dir', help='Maf directory', required=True)
parser.add_argument('-ref', help='Reference species name', required=True)
parser.add_argument('-tree', help='Species tree',
                    default="'((Greattit Groundtit) (Flycatcher Zebrafinch))'")
parser.add_argument('-out', help='Output path and filename', required=True)
args = parser.parse_args()

# variables
maf_dir = args.maf_dir
mafs = maf_dir + '*.maf'
ref_name = args.ref
out_maf = args.out
out_dir = out_maf[:out_maf.rfind('/')+1]
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
temp_dir = out_dir + 'roast_temp'
if not os.path.isdir(temp_dir):
    os.mkdir(temp_dir)
tree = args.tree

# construct and submit command line
roast = ('"roast + '
         'T=' + temp_dir + ' '
         'E=' + ref_name + ' ' +
         tree + ' ' +
         mafs + ' ' +
         out_maf + '"')

qsub = ('python qsub_gen.py '
        '-cmd "cd ' + maf_dir + '" '
        '-cmd ' + roast + ' '
        '-o ' + out_dir + 'roast '
        '-t 168 '
        '-jid roast.sh '
        '-OM q')

subprocess.call(qsub, shell=True)
