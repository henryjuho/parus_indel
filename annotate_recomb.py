#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='/input/directory/vcf', required=True)
parser.add_argument('-ref', help='/path/to/reference/genome.fa', required=True)
parser.add_argument('-poly', help='File containing variable values for polynomial equations for each chromosome',
                    required=True)
parser.add_argument('-out', help='/output/directory/', required=True)
parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
args = parser.parse_args()

# variables
vcf = args.vcf
ref = args.ref
poly = args.poly
out_dir = args.out
evolgen = args.evolgen

# get chromosome list
grep_cmd = 'grep -v ^# ' + vcf + ' | cut -f 1 | uniq'
chromo_list = subprocess.Popen(grep_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0].split('\n')[:-1]
chromo_list = [x for x in chromo_list if x.startswith('chr')]

# loop through chromo list and submit annotation job for each
vcf_outs = []
hold_list = []
for chromo in chromo_list:
    chr_jid = 'recomb_' + chromo + '.sh'
    hold_list.append(chr_jid)
    annotate_vcf_cmd = ('./annotate_recomb_chr.py '
                        '-vcf ' + vcf + ' '
                        '-ref ' + ref + ' '
                        '-poly ' + poly + ' '
                        '-chr ' + chromo + ' '
                        '-out ' + out_dir)
    out_vcf = out_dir + vcf[vcf.rfind('/')+1:].replace('.vcf', '.recomb.' + chromo + '.vcf')
    vcf_outs.append(out_vcf)
    q_sub([annotate_vcf_cmd],
          out=out_dir+'recomb_'+chromo,
          t=8,
          jid=chr_jid,
          evolgen=evolgen)

# submit concat job
cat_cmd = ('./catVCFs.py '
           '-out_vcf ' + out_dir + vcf[vcf.rfind('/')+1:].replace('.vcf', '.recomb.vcf'))
for out_vcf in vcf_outs:
    cat_cmd += ' -vcf ' + out_vcf

q_sub([cat_cmd],
      out=out_dir+'recomb_cat',
      hold=hold_list,
      evolgen=evolgen)
