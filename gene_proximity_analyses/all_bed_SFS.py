#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
import os
from qsub import q_sub


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='vcf file with indels', required=True)
    parser.add_argument('-bed_dir', help='bin beds', required=True)
    parser.add_argument('-call_fa', help='callable sites fasta file', required=True)
    parser.add_argument('-out_dir', help='output directory', required=True)
    parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
    args = parser.parse_args()

    # if out dir doesn't exist, create it
    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir)

    # sfs bed append cmd
    cmd = 'zgrep -w ^{chromo} {bed} | python ~/parus_indel/gene_proximity_analyses/addSFS.py {vcf} {call_fa} > {output}'

    # input bed files
    bed_list = [args.bed_dir + x for x in os.listdir(args.bed_dir) if x.endswith('.bed.gz')]
    for dist_bed in bed_list:
        jid_list = []

        basename = dist_bed.split('/')[-1].replace('.bed.gz', '')
        dist_bed_pysam = pysam.TabixFile(dist_bed)
        out_stem = '{}{}_sfs'.format(args.out_dir, basename)

        for chromo in dist_bed_pysam.contigs:

            output_name = '{}_{}.bed'.format(out_stem, chromo)

            sfs_cmd = cmd.format(chromo=chromo, bed=dist_bed, vcf=args.vcf, call_fa=args.call_fa, output=output_name)
            outs = output_name.replace('.bed', '')
            jid = outs.split('/')[-1] + '.sh'
            q_sub([sfs_cmd], t=72, out=outs, jid=jid, evolgen=args.evolgen)
            jid_list.append(jid)

        # gather cmds
        gather_bed = out_stem + '.bed.gz'
        gather = 'cat {}*.bed | sort -k1,1 -k2,2n | bgzip -c > {}'.format(out_stem, gather_bed)
        tabix = 'tabix -pbed {}'.format(gather_bed)
        q_sub([gather, tabix], out=out_stem, hold=jid_list, evolgen=args.evolgen)


if __name__ == '__main__':
    main()
