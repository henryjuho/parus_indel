#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
from vcf2raw_sfs import vcf2sfs
import anavar_utils as an
from qsub import q_sub
import gzip
import random
import sys
import os
import subprocess
sys.path.insert(0, os.getenv('HOME') + '/parus_indel/anavar_analyses')
from sel_vs_neu_anavar import sfs2counts, read_callable_csv
from window_sel_vs_neu_anavar import window_call_sites


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='vcf file with indels', required=True)
    parser.add_argument('-bed_list', help='list of bin beds', required=True)
    parser.add_argument('-call_fa', help='callable sites fasta file', required=True)
    parser.add_argument('-call_csv', default='/fastdata/bop15hjb/GT_ref/gt_callable_summary.csv',
                        help=argparse.SUPPRESS)
    parser.add_argument('-equal_theta', help='if specified runs with equal mutation ratesbetween neu and sel sites',
                        default=False, action='store_true')
    parser.add_argument('-dfe', help='determines type of distribution to have in model',
                        choices=['discrete', 'continuous'], default='continuous')
    parser.add_argument('-out_pre', help='output prefix', required=True)
    parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
    parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
    args = parser.parse_args()

    # submission loop
    if args.sub is True:
        command_line = [' '.join([x for x in sys.argv if x != '-sub'])]
        q_sub(command_line, args.out_pre + '.control_job')
        sys.exit()

    # flags
    if args.equal_theta:
        constraint = 'equal_mutation_rate'
    else:
        constraint = 'none'

    # files
    call_fasta = pysam.FastaFile(args.call_fa)

    # make a sorted list of form [('1', 'bin1.bed.gz'), ('2', 'bin2.bed.gz')]
    bed_files = sorted([(x.rstrip().split('.')[-3].replace('bin', ''), x.rstrip()) for
                        x in open(args.bed_list) if x.rstrip().endswith('.bed.gz')])

    # neu reference
    n_i_sfs = vcf2sfs(vcf_name=args.vcf, mode='ins',
                      auto_only=True,
                      regions=['intergenic_ar', 'intron_ar'])

    n_d_sfs = vcf2sfs(vcf_name=args.vcf, mode='del',
                      auto_only=True,
                      regions=['intergenic_ar', 'intron_ar'])

    sfs_ni = sfs2counts(n_i_sfs, 20)
    sfs_nd = sfs2counts(n_d_sfs, 20)

    neu_m = read_callable_csv(args.call_csv)['ALL']['AR']['pol']

    # everything else per window
    for dist_bin in bed_files:

        bin_id = dist_bin[0]
        bin_bed = dist_bin[1]

        out_stem = '{}_bin{}'.format(args.out_pre, bin_id)

        call_sites = 0

        sfs_cmd = ('bedtools intersect -header -a {} -b {} | '
                   '~/sfs_utils/vcf2raw_sfs.py -region intergenic -region intron -mode {} -auto_only')

        ins_sfs = subprocess.Popen(sfs_cmd.format(args.vcf, bin_bed, 'ins'), shell=True, stdout=subprocess.PIPE
                                   ).communicate()[0].split('\n')[:-1]
        del_sfs = subprocess.Popen(sfs_cmd.format(args.vcf, bin_bed, 'del'), shell=True, stdout=subprocess.PIPE
                                   ).communicate()[0].split('\n')[:-1]

        for coord_set in gzip.open(bin_bed):

            coords = coord_set.split()

            if coords[0] == 'chrZ':
                continue

            # get sel call sites
            call_sites += window_call_sites(call_fasta, None, (coords[0], int(coords[1]), int(coords[2])))

        ins_counts = sfs2counts(ins_sfs, 20)
        del_counts = sfs2counts(del_sfs, 20)

        # anavar setup
        sfs_data = {'selected_INS': (ins_counts, call_sites), 'selected_DEL': (del_counts, call_sites),
                    'neutral_INS': (sfs_ni, neu_m), 'neutral_DEL': (sfs_nd, neu_m)}

        anavar_path = '/shared/evolgen1/shared_data/program_files/sharc/'

        anavar_cmd = '{path}anavar1.22 {ctl} {rslts} {log} {seed}'

        # sort file names
        ctl_name = out_stem + '.control.txt'
        result_name = out_stem + '.results.txt'
        log_name = out_stem + '.log.txt'

        # make control file
        ctl = an.IndelNeuSelControlFile()

        ctl.set_alg_opts(search=500, alg='NLOPT_LD_SLSQP', key=3,
                         epsabs=1e-20, epsrel=1e-9, rftol=1e-9,
                         maxtime=3600, optional=True)

        ctl.set_data(sfs_data, 20, dfe=args.dfe, c=1,
                     theta_r=(1e-10, 0.1), r_r=(0.01, 100),
                     scale_r=(0.1, 5000.0), gamma_r=(-5e4, 1e2))

        ctl.set_constraint(constraint)
        ctl_contents = ctl.construct()

        with open(ctl_name, 'w') as control:
            control.write(ctl_contents)

        # submit anavar window job
        window_cmd = anavar_cmd.format(path=anavar_path, ctl=ctl_name, rslts=result_name, log=log_name,
                                       seed=random.randint(0, 1000000))

        q_sub([window_cmd], out=out_stem, t=24, evolgen=args.evolgen)


if __name__ == '__main__':
    main()
