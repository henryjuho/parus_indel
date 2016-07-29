#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gff', help='GFF annotation file', required=True)
parser.add_argument('-vcf', help='VCF file to annotate variants in', required=True)
parser.add_argument('-ref', help='Reference genome', required=True)
parser.add_argument('-db_dir', help='Location of databases', required=True)
parser.add_argument('-out', help='Output directory', required=True)
parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
args = parser.parse_args()

# variables
gff = args.gff
vcf = args.vcf
ref = args.ref
db_dir = args.db_dir
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
    chr_jid = 'degeneracy_annotation_' + chromo + '.sh'
    hold_list.append(chr_jid)
    out_vcf = out_dir + vcf[vcf.rfind('/')+1:].replace('.vcf', '.degen.' + chromo + '.vcf')

    # command to annotate chromosomal vcf
    annotate_degen_cmd = ('./annotate_degen_chr_2.py '
                          '-gff ' + gff + ' '
                          '-vcf ' + vcf + ' '
                          '-ref ' + ref + ' '
                          '-chr ' + chromo + ' '
                          '-db_dir ' + db_dir + ' '
                          '-out ' + out_dir)

    vcf_outs.append(out_vcf)
    q_sub([annotate_degen_cmd],
          out=out_dir + 'degeneracy_annotation_' + chromo,
          t=48,
          jid=chr_jid,
          evolgen=evolgen)

# submit concat job
cat_cmd = ('./catVCFs.py '
           '-out_vcf ' + out_dir + vcf[vcf.rfind('/')+1:].replace('.vcf', '.degen.vcf'))
for out_vcf in vcf_outs:
    cat_cmd += ' -vcf ' + out_vcf

q_sub([cat_cmd],
      out=out_dir + 'vcf_concatenation',
      hold=hold_list,
      evolgen=evolgen)
