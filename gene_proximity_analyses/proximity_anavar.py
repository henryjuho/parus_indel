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
from itertools import chain
sys.path.insert(0, os.getenv('HOME') + '/parus_indel/anavar_analyses')
from sel_vs_neu_anavar import sfs2counts, read_callable_csv
from window_sel_vs_neu_anavar import window_call_sites


def bootstrap(boot_list_list):

    new_list_list = [[] for y in boot_list_list]

    for i in range(0, len(boot_list_list[0])):
        i = random.randint(0, len(boot_list_list[0])-1)

        for x in range(0, len(boot_list_list)):

            new_list_list[x].append(boot_list_list[x][i])

    return new_list_list


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
    parser.add_argument('-bootstrap', help='Number of rounds of bootstrapping to perform', default=0, type=int)
    parser.add_argument('-out_pre', help='output prefix', required=True)
    parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
    parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
    args = parser.parse_args()

    # submission loop
    if args.sub is True:
        command_line = [' '.join([x for x in sys.argv if x != '-sub'])]
        q_sub(command_line, args.out_pre + '.control_job', evolgen=args.evolgen)
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

    anavar_path = '/shared/evolgen1/shared_data/program_files/sharc/'

    anavar_cmd = '{path}anavar1.22 {ctl} {rslts} {log} {seed}'

    # everything else per window
    seed = 0
    for dist_bin in bed_files:
        seed += 1

        bin_id = dist_bin[0]
        bin_bed = dist_bin[1]

        out_stem = '{}_bin{}'.format(args.out_pre, bin_id)

        bootstrapable_data = {'del': [], 'ins': [], 'call': []}

        sfs_cmd = ('tabix -h {vcf} {chromo}:{start}-{end} | '
                   '~/sfs_utils/vcf2raw_sfs.py -mode {mode}')

        for coord_set in gzip.open(bin_bed):

            coords = coord_set.split()

            if coords[0] == 'chrZ':
                continue

            for indel_type in ['ins', 'del']:

                sfs = subprocess.Popen(sfs_cmd.format(
                    chromo=coords[0], start=coords[1], end=coords[2],
                    vcf=args.vcf, mode=indel_type), shell=True, stdout=subprocess.PIPE
                                   ).communicate()[0].split('\n')[:-1]

                bootstrapable_data[indel_type].append(sfs)

            # get sel call sites
            bootstrapable_data['call'].append(window_call_sites(call_fasta, None,
                                                                (coords[0], int(coords[1]), int(coords[2]))))

        # ====BOOTSTRAPPING=====
        # sort file names
        sge_tag = '$SGE_TASK_ID'
        ctl_template = '{}_bootstrap{}.control.txt'

        ctl_name = ctl_template.format(out_stem, sge_tag)
        result_name = '{}_bootstrap$SGE_TASK_ID.results.txt'.format(out_stem)
        log_name = '{}_bootstrap$SGE_TASK_ID.log.txt'.format(out_stem)

        for i in range(1, args.bootstrap+2):

            ins_freqs = bootstrapable_data['ins']
            del_freqs = bootstrapable_data['del']
            call_list = bootstrapable_data['call']

            if i == 1:
                ins_sfs = list(chain.from_iterable(ins_freqs))
                del_sfs = list(chain.from_iterable(del_freqs))
                call_sites = sum(call_list)

            else:

                bs_ins, bs_del, bs_call = bootstrap([ins_freqs, del_freqs, call_list])

                ins_sfs = list(chain.from_iterable(bs_ins))
                del_sfs = list(chain.from_iterable(bs_del))
                call_sites = sum(bs_call)

            ins_counts = sfs2counts(ins_sfs, 20)
            del_counts = sfs2counts(del_sfs, 20)

            # anavar setup
            sfs_data = {'selected_INS': (ins_counts, call_sites), 'selected_DEL': (del_counts, call_sites),
                        'neutral_INS': (sfs_ni, neu_m), 'neutral_DEL': (sfs_nd, neu_m)}

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

            with open(ctl_template.format(out_stem, i), 'w') as control:
                control.write(ctl_contents)

        # submit anavar window job
        window_cmd = anavar_cmd.format(path=anavar_path, ctl=ctl_name, rslts=result_name, log=log_name,
                                       seed=seed)

        q_sub([window_cmd], out=out_stem, t=24, evolgen=args.evolgen, array=[1, args.bootstrap + 1])


if __name__ == '__main__':
    main()
