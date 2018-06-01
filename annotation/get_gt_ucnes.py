#!/usr/bin/env python

from __future__ import print_function
import argparse
import subprocess
from qsub import q_print as q_sub


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-wga', help='wga bed file', required=True)
    parser.add_argument('-ucne_bed', help='UCNE bed file', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
    args = parser.parse_args()

    # get chromosome list
    grep_cmd = 'zcat {} | cut -f 1 | uniq'.format(args.ucne_bed)
    chromo_list = subprocess.Popen(grep_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0].split('\n')[:-1]
    chromo_list = [x for x in chromo_list if x.startswith('chr')]

    jids = []

    for chromo in chromo_list:

        out = 'gt_ucne_{}.bed'.format(chromo)

        cmd = ('zcat {} | '
               '~/WGAbed/non_ref_intersect.py '
               '-b {} -q Zebrafinch -c {} | '
               'cut -f 1-3 | '
               'sort -k1,1 -k2,2n | '
               'bedtools merge '
               '> {}').format(args.wga, args.ucne_bed, chromo, args.out_dir + out)

        jid = out.replace('.bed', '.sh')
        jids.append(jid)

        q_sub([cmd], out=args.out_dir + out.replace('.bed', ''), jid=jid, evolgen=args.evolgen)

    # gather
    gather = 'zcat {}*.bed | bgzip -c > {}gt_ucne.bed.gz'.format(args.out_dir, args.out_dir)
    index = 'tabix -pbed {}gt_ucne.bed.gz'.format(args.out_dir)

    q_sub([gather, index], out=args.out_dir + 'ucne_bed_merge', hold=jids, evolgen=args.evolgen)


if __name__ == '__main__':
    main()
