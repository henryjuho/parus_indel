#!/usr/bin/env python

import argparse
import subprocess
from hen_utils import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf',
                    help='VCF file to get site frequency spectrum from',
                    required=True)
parser.add_argument('-folded',
                    help='Specify Y for folded and N for unfolded spectrum',
                    default='Y',
                    choices=['Y', 'N'])
parser.add_argument('-bin',
                    help='Specify genomic regions to bin data by (requires an annotated vcf',
                    default=[],
                    choices=['CDS_non_frameshift', 'CDS_frameshift', 'intron', 'intergenic', 'CDS',
                             'intergenic_ww_ss', 'zerofold'],
                    action='append')
parser.add_argument('-mode', help='Mode to run, SNP or INDEL', choices=['SNP', 'INDEL'], required=True)
parser.add_argument('-sfs_out',
                    help='Output sfs data prefix',
                    required=True)
parser.add_argument('-bootstrap', help='Help specify number of times to bootstrap, default is none', default='0')
parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
args = parser.parse_args()

# variables
vcf_file = args.vcf
folded = args.folded
output = args.sfs_out
outdir = file_location(output)
bins = args.bin
bootstrap = args.bootstrap
mode = args.mode
evolgen = args.evolgen

if mode == 'SNP':
    sfs_script = './snpSFS.py '
else:
    sfs_script = './indelSFS.py '

# get chromosomal vcfs
chromosomes = subprocess.Popen('zgrep -v ^# ' + vcf_file + ' | cut -f 1 | uniq', shell=True,
                               stdout=subprocess.PIPE).communicate()[0].split()

for chromo in chromosomes:

    # generate chromosomal vcf
    chromo_vcf = outdir + extract_filename(vcf_file).rstrip('.gz').rstrip('.vcf') + '.' + chromo + '.vcf'
    subprocess.call('zgrep ^# ' + vcf_file + ' > ' + chromo_vcf, shell=True)
    subprocess.call('zgrep ^' + chromo + ' ' + vcf_file + ' >> ' + chromo_vcf, shell=True)

    # output
    chromo_output = output + '.' + chromo

    # run SFS script
    sfs_cmd = (sfs_script +
               '-vcf ' + chromo_vcf + ' '
               '-folded ' + folded + ' ')

    sfs_cmd += ''.join(['-bin ' + b + ' ' for b in bins])

    sfs_cmd += ('-include_sex '
                '-bootstrap ' + bootstrap + ' '
                '-sfs_out ' + chromo_output + ' '
                '-sub ')

    if evolgen is True:
        sfs_cmd += '-evolgen'

    subprocess.call(sfs_cmd, shell=True)
