#!/usr/bin/env python

import argparse
import os
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-dir', help='Directory containing maf files to ensure single coverage for', required=True)
args = parser.parse_args()

# variables
directory = args.dir
maf_list = [maf for maf in os.listdir(directory) if maf.endswith('.maf')]
out_dir = directory + 'single_coverage/'
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

# run single_cov2
for maf in maf_list:
    in_maf = directory + maf
    output = out_dir + maf.rstrip('.maf') + '.sing.maf'
    cmd_line = ('"single_cov2 ' +
                in_maf + ' [R=Chicken] > ' +
                output + '"')
    qsub_cmd = ('python qsub_gen.py '
                '-cmd ' + cmd_line + ' '
                '-o ' + out_dir + 'single_cov2_' + maf.rstrip('.maf') + ' '
                '-jid single_cov2_' + maf.rstrip('.maf') + '.sh '
                '-OM q')
    subprocess.call(qsub_cmd, shell=True)
