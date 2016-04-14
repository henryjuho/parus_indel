#!/usr/bin/env python

import argparse
import subprocess
import os

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-ref', help='Whole genome fasta of reference species', required=True)
parser.add_argument('-ref_name', help='Name of ref spp for file naming', required=True)
parser.add_argument('-fa_list', help='List of chromosomal fastas to align to ref', required=True)
parser.add_argument('-out', help='Output directory', required=True)
args = parser.parse_args()

# variables
query_genomes = [fa.rstrip('\n').split() for fa in open(args.fa_list)]
ref_genome = args.ref
ref_name = args.ref_name
out_dir = args.out
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)

# iterate through query genomes
for genome in query_genomes:
    output_prefix = out_dir + ref_name + '.' + genome[1]

    # construct commandline
    lastz_cmd = ('"lastz_latest ' +
                 ref_genome + '[multiple][nameparse=darkspace] ' +
                 genome[0] + '[nickname=' + genome[1] + '] '
                 '--gfextend '
                 '--chain '
                 '--gapped '
                 '--rdotplot=' + output_prefix + '.rdotplot '
                 '--format=maf '
                 '> ' + output_prefix + '.maf"')

    r_plotting_cmd = ('"Rscript align_dots.R ' +
                      output_prefix + '.rdotplot ' +
                      output_prefix + '.dotplot.pdf"')

    # submit to SGE
    qsub_cmd = ('python qsub_gen.py '
                '-cmd ' + lastz_cmd + ' '
                '-cmd ' + r_plotting_cmd + ' '
                '-t 72 '
                '-mem 10 -rmem 10 '
                '-o ' + output_prefix + ' '
                '-OM q')
    subprocess.call(qsub_cmd, shell=True)
