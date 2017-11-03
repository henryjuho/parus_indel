#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
from vcf2raw_sfs import vcf2sfs
from sel_vs_neu_anavar import sfs2counts, read_callable_csv
import anavar_utils as an
from qsub import q_sub


def window_call_sites(call_fa, region_bed, window_coords):

    """
    returns number of callable sites for specied region in window
    :param call_fa: pysam.FastaFile()
    :param region_bed: pysam.TabixFile()
    :param window_coords: tuple
    :return: int
    """

    regions = region_bed.fetch(window_coords[0], window_coords[1], window_coords[2], parser=pysam.asTuple())

    call_sites = 0

    for reg in regions:
        call_seq = call_fa.fetch(reg[0], int(reg[1]), int(reg[2]))
        call_sites += call_seq.count('K')

    return call_sites


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='vcf file with indels', required=True)
    parser.add_argument('-windows', help='window file', required=True)
    parser.add_argument('-call_fa', help='callable sites fasta file', required=True)
    parser.add_argument('-noncoding_bed', help='bed file of non-coding regions', required=True)
    parser.add_argument('-call_csv', default='/fastdata/bop15hjb/GT_ref/gt_callable_summary.csv',
                        help=argparse.SUPPRESS)
    parser.add_argument('-equal_theta', help='if specified runs with equal mutation ratesbetween neu and sel sites',
                        default=False, action='store_true')
    parser.add_argument('-out_pre', help='output prefix', required=True)
    parser.add_argument('-evolgen', help='if specified will run on lab queue', default=False, action='store_true')
    args = parser.parse_args()

    # flags
    if args.equal_theta:
        constraint = 'equal_mutation_rate'
    else:
        constraint = 'none'

    # files
    call_fasta = pysam.FastaFile(args.call_fa)
    nc_bed = pysam.TabixFile(args.noncoding_bed)

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
    for window in open(args.windows):

        window_info = window.split()
        coords = window_info[0:3]
        window_id = window_info[3]
        out_stem = '{}_window{}'.format(args.out_pre, window_id)

        # get sel sfs
        ins_sfs = vcf2sfs(args.vcf, mode='ins', regions=['intron', 'intergenic'],
                          chromo=coords[0], start=int(coords[1]), stop=int(coords[2]))
        del_sfs = vcf2sfs(args.vcf, mode='del', regions=['intron', 'intergenic'],
                          chromo=coords[0], start=int(coords[1]), stop=int(coords[2]))

        ins_counts = sfs2counts(ins_sfs, 20)
        del_counts = sfs2counts(del_sfs, 20)

        # get sel call sites
        call_sites = window_call_sites(call_fasta, nc_bed, coords)

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

        ctl.set_data(sfs_data, 20, dfe='continuous', c=1,
                     theta_r=(1e-10, 0.1), r_r=(0.01, 100), scale_r=(0.1, 5000.0))

        ctl.set_constraint(constraint)
        ctl_contents = ctl.construct()

        with open(ctl_name, 'w') as control:
            control.write(ctl_contents)

        # submit anavar window job
        window_cmd = anavar_cmd.format(path=anavar_path, ctl=ctl_name, rslts=result_name, log=log_name,
                                       seed=int(window_id))

        q_sub([window_cmd], out=out_stem, jid=out_stem + '.sh', t=8, evolgen=args.evolgen)


if __name__ == '__main__':
    main()
