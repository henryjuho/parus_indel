#!/usr/bin/env python

from __future__ import print_function
import anavar_utils as an
import argparse
from qsub import q_sub
import random
from collections import Counter
from vcf2raw_sfs import vcf2sfs


def read_callable_csv(csv):

    """
    reads in a callable sites summary file
    :param csv: str
    :return: dict
    """

    call_sites = {}
    for line in open(csv):
        if not line.startswith('contig'):
            info = line.rstrip().split(',')
            contig, reg, all_call, pol_call = info
            if contig not in call_sites.keys():
                call_sites[contig] = {reg: {'all': float(all_call), 'pol': float(pol_call)}}
            else:
                call_sites[contig][reg] = {'all': float(all_call), 'pol': float(pol_call)}
    return call_sites


def read_callable_txt(txt):

    """
    reads in a callable sites summary file
    :param txt: str
    :return: dict
    """

    call_sites = {'ALL': {}}
    for line in open(txt):
        if not line.startswith('category'):
            info = line.rstrip().split('\t')
            reg, var, call = info

            if reg not in call_sites['ALL'].keys():
                call_sites['ALL'][reg] = {}

            if var == 'SNP':
                call_sites['ALL'][reg]['all'] = float(call)

            elif var == 'INS':
                call_sites['ALL'][reg]['pol'] = float(call)

            else:
                continue

    return call_sites


def resample_replace(site_freqs):

    """
    resamples sfs with replacement
    :param site_freqs: list
    :return: list
    """

    resamp_sfs = []
    for i in range(0, len(site_freqs)):
        random_no = random.randint(0, len(site_freqs) - 1)
        resamp_sfs.append(site_freqs[random_no])
    return resamp_sfs


def sfs2counts(freq_list, n):

    """
    converts list of sfs to a condensed sfs that can be plotted or passed to anavar etc
    :param freq_list: list
    :param n: int
    :return: list
    """

    pos_biallelic_freqs = [round(i / float(n), 3) for i in range(1, int(n))]

    counted_sorted_sfs = sorted(Counter([str(x) for x in freq_list]).most_common(), key=lambda z: z[0])
    sfs_freq_dict = {x[0]: x[1] for x in counted_sorted_sfs}

    counts = []
    for frequency in pos_biallelic_freqs:
        try:
            no_var = sfs_freq_dict[str(frequency)]
        except KeyError:
            no_var = 0  # account for freqs with 0 variants

        counts.append(no_var)

    return counts


def prepare_indel_sfs(vcf, call, n, sel_sfs_regions, call_sel_reg):

    """
    gets sfs from vcf and prepares as anavar input
    :param vcf: str
    :param call: dict
    :param n: int
    :param sel_sfs_regions: list
    :param call_sel_reg: str
    :return: dict
    """

    # extract site frequencies
    del_sfs = vcf2sfs(vcf_name=vcf, mode='del',
                      auto_only=True,
                      regions=sel_sfs_regions)

    n_d_sfs = vcf2sfs(vcf_name=vcf, mode='del',
                      auto_only=True,
                      regions=['intergenic_ar', 'intron_ar'])

    ins_sfs = vcf2sfs(vcf_name=vcf, mode='ins',
                      auto_only=True,
                      regions=sel_sfs_regions)

    n_i_sfs = vcf2sfs(vcf_name=vcf, mode='ins',
                      auto_only=True,
                      regions=['intergenic_ar', 'intron_ar'])

    # convert to correct format for anavar
    sfs_i = sfs2counts(ins_sfs, n)
    sfs_d = sfs2counts(del_sfs, n)
    sfs_ni = sfs2counts(n_i_sfs, n)
    sfs_nd = sfs2counts(n_d_sfs, n)

    # get callable sites
    sel_m = call['ALL'][call_sel_reg]['pol']
    neu_m = call['ALL']['AR']['pol']

    # construct control file sfs
    sfs_m = {'selected_INS': (sfs_i, sel_m), 'selected_DEL': (sfs_d, sel_m),
             'neutral_INS': (sfs_ni, neu_m), 'neutral_DEL': (sfs_nd, neu_m)}

    return sfs_m


def prepare_snp_sfs(vcf, call, n, sel_sfs_regions, call_sel_reg):

    """
    gets sfs from vcf and prepares as anavar input
    :param vcf: str
    :param call: dict
    :param n: int
    :param sel_sfs_regions: list
    :param call_sel_reg: str
    :return: dict
    """

    # extract site frequencies
    sel_sfs = vcf2sfs(vcf_name=vcf, mode='snp',
                      auto_only=True,
                      regions=sel_sfs_regions)

    neu_sfs = vcf2sfs(vcf_name=vcf, mode='snp',
                      auto_only=True,
                      degen=4)

    # convert to correct format for anavar
    sfs_sel = sfs2counts(sel_sfs, n)
    sfs_neu = sfs2counts(neu_sfs, n)

    # get callable sites
    sel_m = call['ALL'][call_sel_reg]['pol']
    neu_m = call['ALL']['4fold']['pol']

    # construct control file sfs
    sfs_m = {'selected_SNP': (sfs_sel, sel_m), 'neutral_SNP': (sfs_neu, neu_m)}

    return sfs_m


def sel_v_neu_anavar(mode, vcf, call, sel_region, constraint, n, c, dfe, alg, nnoimp, maximp,
                     out_stem, search, degree, spread, evolgen, start_index):

    """
    submits anavar jobs to cluster after writing required files etc
    :param mode: str
    :param vcf: str
    :param call: dict
    :param sel_region: str
    :param constraint: str
    :param n: int
    :param c: int
    :param dfe: str
    :param alg: str
    :param nnoimp: int
    :param maximp: int
    :param out_stem: str
    :param search: int
    :param degree: int
    :param spread: int
    :param evolgen: bool
    :param start_index: int
    :return: None
    """

    anavar_path = '/shared/evolgen1/shared_data/program_files/sharc/'

    anavar_cmd = '{path}anavar1.22 {ctl} {rslts} {log} {seed}'

    # sort file names
    ctl_name = out_stem + '.control.txt'

    # region combinations
    region_combs = {'CDS': ['CDS_frameshift', 'CDS_non_frameshift'],
                    'intron': ['intron'],
                    'intergenic': ['intergenic'],
                    'noncoding': ['intergenic', 'intron']}

    # make control file
    if mode == 'snp':
        sfs_data = prepare_snp_sfs(vcf, call, n,
                                   sel_sfs_regions=region_combs[sel_region],
                                   call_sel_reg=sel_region)
        ctl = an.SNPNeuSelControlFile()

    else:
        sfs_data = prepare_indel_sfs(vcf, call, n,
                                     sel_sfs_regions=region_combs[sel_region],
                                     call_sel_reg=sel_region)
        ctl = an.IndelNeuSelControlFile()

    ctl.set_alg_opts(search=search, alg=alg, key=3,
                     epsabs=1e-20, epsrel=1e-9, rftol=1e-9,
                     maxtime=3600, optional=True,
                     maximp=maximp, nnoimp=nnoimp)

    ctl.set_data(sfs_data, n, dfe=dfe, c=c, gamma_r=(-5e4, 1e5), theta_r=(1e-12, 0.1), r_r=(0.01, 100),
                 scale_r=(0.1, 5000.0))
    if degree != 50:
        ctl.set_dfe_optional_opts(degree=degree, optional=True)
    ctl.set_constraint(constraint)
    ctl_contents = ctl.construct()
    with open(ctl_name, 'w') as control:
        control.write(ctl_contents)

    res_file_list = out_stem + '.allres.list.txt'
    hjids = []
    with open(res_file_list, 'a') as res_list:

        # split into requested jobs
        for i in range(start_index, start_index+spread):

            split_stem = '{}.split{}'.format(out_stem, i)

            result_name = split_stem + '.results.txt'
            log_name = split_stem + '.log.txt'

            print(result_name, file=res_list)

            # call anavar
            rep_cmd = anavar_cmd.format(path=anavar_path, ctl=ctl_name, rslts=result_name, log=log_name, seed=i)

            q_sub([rep_cmd], out=split_stem, jid=split_stem.split('/')[-1] + '.sh', t=8, evolgen=evolgen)
            hjids.append(split_stem.split('/')[-1] + '.sh')

    # hold job to merge outputs
    merge_out = out_stem + '.merged.results.txt'
    gather = 'cat {} | gather_searches.py {}'.format(res_file_list, merge_out)
    q_sub([gather], out=out_stem + '.merge', hold=hjids, evolgen=evolgen)


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-mode', help='variant type to run on', choices=['snp', 'indel'], required=True)
    parser.add_argument('-call_txt', help='Callable sites summary file', required=True)
    parser.add_argument('-vcf', help='VCF file to extract site frequencies from', required=True)
    parser.add_argument('-n', help='Sample size', required=True)
    parser.add_argument('-sel_type', help='Subset of variants to have as sel sfs',
                        choices=['CDS', 'intron', 'intergenic', 'noncoding'],
                        default='CDS')
    parser.add_argument('-c', help='Number of classes to run model with', required=True, type=int)
    parser.add_argument('-dfe', help='type of dfe to fit, discrete or continuous', default='discrete',
                        choices=['discrete', 'continuous'])
    parser.add_argument('-constraint', help='Constraint for model', choices=['none', 'equal_mutation_rate'],
                        default='none')
    parser.add_argument('-n_search', help='Number of searches to conduct per job', default=500, type=int)
    parser.add_argument('-alg', help='Algorithm to use', default='NLOPT_LD_SLSQP',
                        choices=['NLOPT_LD_SLSQP', 'NLOPT_LD_LBFGS', 'NLOPT_LN_NELDERMEAD'])
    parser.add_argument('-nnoimp', help='nnoimp value', default=1, type=int)
    parser.add_argument('-maximp', help='maximp value', default=3, type=int)
    parser.add_argument('-split', help='Number of jobs to split runs across, each job will run the control file once'
                                       'with a different seed given to anavar', default=1, type=int)
    parser.add_argument('-degree', help='changes degree setting in anavar', default=50, type=int)
    parser.add_argument('-out_pre', help='File path and prefix for output', required=True)
    parser.add_argument('-start_index', help='ID for first bin and for first seed', default=1, type=int)
    parser.add_argument('-evolgen', help='If specified will run on evolgen', default=False, action='store_true')
    args = parser.parse_args()

    # variables
    call_site_dict = read_callable_txt(args.call_txt)
    out_pre = args.out_pre

    # construct process
    sel_v_neu_anavar(mode=args.mode,
                     vcf=args.vcf, call=call_site_dict,
                     sel_region=args.sel_type,
                     constraint=args.constraint,
                     n=args.n, c=args.c, dfe=args.dfe,
                     alg=args.alg,
                     nnoimp=args.nnoimp, maximp=args.maximp,
                     out_stem=out_pre,
                     search=args.n_search,
                     degree=args.degree,
                     spread=args.split,
                     evolgen=args.evolgen,
                     start_index=args.start_index)

if __name__ == '__main__':
    main()
