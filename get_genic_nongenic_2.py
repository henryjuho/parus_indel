#!/usr/bin/env python

import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF file with variants you want to subset', required=True)
parser.add_argument('-bed', help='Bed file with annotation data in', required=True)
parser.add_argument('-out', help='Output directory', required=True)
args = parser.parse_args()

# input files
vcf_file = args.vcf
bed_file = args.bed

# output locations
out_dir = args.out
gene_vcf = out_dir + vcf_file[vcf_file.rfind('/')+1:].rstrip('.vcf') + '.genic.vcf'
non_genic_vcf = out_dir + vcf_file[vcf_file.rfind('/')+1:].rstrip('.vcf') + '.intergenic.vcf'

# construct bedtools command lines
genic = ('"bedtools intersect '
         '-header '
         '-wa '
         '-a ' + vcf_file + ' '
         '-b ' + bed_file + ' '
         '> ' + gene_vcf + '"')

non_genic = ('"bedtools intersect '
             '-header '
             '-wa '
             '-v '
             '-a ' + vcf_file + ' '
             '-b ' + bed_file + ' '
             '> ' + non_genic_vcf + '"')

qsubber = ('python qsub_gen.py '
           '-cmd ' + genic + ' '
           '-cmd ' + non_genic + ' '
           '-o ' + out_dir + 'get_genic_regions '
           '-OM q '
           '-evolgen')

subprocess.call(qsubber, shell=True)
