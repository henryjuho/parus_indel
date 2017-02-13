#!/usr/bin/env python

import argparse
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-wga_bed', help='wga.bed alignment file', required=True)
parser.add_argument('-r_bed', help='bed file with LINE coordinates for reference species', required=True)
parser.add_argument('-q_bed', help='spp name and bed file with LINE coordinates, comma separated: spp,file '
                                   'specify for each non_ref species in alignment',
                    action='append', required=True)
parser.add_argument('-chr_list', help='text file listing chromosomes to extract', required=True)
parser.add_argument('-out', help='output directory and file prefix', required=True)
parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
args = parser.parse_args()

# variables
wga = args.wga_bed
ref_lines = args.r_bed
query_lines = [y.split(',') for y in args.q_bed]
chromo_list = [x.split()[0] for x in open(args.chr_list)]
out = args.out
evolgen = args.evolgen

# chromo loop
out_list = []
jid_list = []
for chromo in chromo_list:

    # intersect commands
    chromo_out = out + '.LINEs.' + chromo + '.wga.bed.gz'
    out_list.append(chromo_out)
    threeway_intersect = ('zgrep ^' + chromo + ' ' + wga + ' | '
                          'bedtools intersect -a stdin -b ' + ref_lines + ' ')
    for query in query_lines:
        threeway_intersect += '| non_ref_intersect.py -b ' + query[1] + ' -q ' + query[0] + ' -c ' + chromo + ' '
    threeway_intersect += ' | bgzip -c > ' + chromo_out

    # submit to cluster
    jid = chromo + '.LINEs.extract.sh'
    jid_list.append(jid)
    q_sub([threeway_intersect], out=out + '.' + chromo + '.LINEs', jid=jid, evolgen=evolgen)

# concat
cat = 'zcat '
for bed in out_list:
    cat += bed + ' '
cat += ' | bgzip -c > ' + out + '.LINEs.wga.bed.gz'
q_sub([cat], out=out + '.conservedLINEs', hold=jid_list, evolgen=evolgen)
