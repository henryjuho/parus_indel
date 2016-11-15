#!/usr/bin/env python

import argparse
from qsub import *
import math
import pysam
import gffutils
import re
import sys
import random
import numpy

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='Vcf file to get summary stats for', required=True)
parser.add_argument('-call_fa', help='Fasta file of callable sites, coded as 0=N, 1=notcalled and 2=called',
                    required=True)
parser.add_argument('-mode', help='Variant mode, either SNP, SNP_comp or INDEL', required=True,
                    choices=['SNP', 'SNP_comp', 'INDEL'])
parser.add_argument('-by_region', help='Specify location of gff databases for use with gffutils', default='NONE')
parser.add_argument('-bootstrap', help='If specified will perform the requested number opf rounds of bootstrapping',
                    default=0)
parser.add_argument('-md', help='If specified will output in markdown table format', default=False, action='store_true')
parser.add_argument('-out', help='Output directory and file name, only required in conjunction with sub', default='')
parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
parser.add_argument('-evolgen', help='If specified will run on lab queue', action='store_true', default=False)
args = parser.parse_args()


# variables
vcf_file = args.vcf
call_fasta = pysam.FastaFile(args.call_fa)
mode = args.mode
markdown = args.md
database_path = args.by_region
bootstrap = int(args.bootstrap)
chromo_list = []

if mode == 'SNP_comp' and database_path == 'NONE':
    sys.exit('-mode SNP_comp can only be run in conjunction with -by_region')

if len(args.out) == 0 and args.sub is True:
    sys.exit('-sub must be specified in conjunction with -out')

# submission loop
if args.sub is True:
    command_line = [' '.join([x for x in sys.argv if x != '-sub' and x != '-evolgen']) + ' > ' + args.out]
    q_sub(command_line, args.out.rstrip('.txt'), evolgen=args.evolgen)
    sys.exit()

# preset dictionaries
if database_path == 'NONE':
    regions = ['NONE']
    database_dict = {}
else:
    regions = ['CDS', 'intergenic', 'intron']
    database_dict = {q.split('.')[-2]: args.by_region + q for q in os.listdir(args.by_region) if q.endswith('.db')}


# functions
def theta_w(n, segsites):
    # theta = S/a
    alpha = sum(1.0/z for z in range(1, n))
    theta = float(segsites) / alpha
    return theta


def pi(n, allele_frq_list):
    no_seqs_dif = [(1.0 - raf**2 - (1.0-raf)**2) * (n/(n-1.0)) for raf in allele_frq_list]
    seq_pi = sum(no_seqs_dif)
    return seq_pi


def tajimas_d(n, allele_frq_list):
    segsites = float(len(allele_frq_list))
    little_d = pi(n, allele_frq_list) - theta_w(n, segsites)

    a1 = sum(1.0/z for z in range(1, n))
    a2 = sum(1.0/z**2 for z in range(1, n))

    e1 = (1.0 / a1) * (((n + 1.0) / (3.0 * (n - 1.0))) - (1.0 / a1))
    e2 = (1.0 / (a1**2 + a2)) * \
         (((2.0 * (n**2 + n + 3.0)) / ((9.0 * n) * (n - 1.0))) -
          ((n + 2.0) / (n * a1)) +
          (a2 / a1**2))
    # print(e2)
    vd = (e1 * segsites) + ((e2 * segsites) * (segsites - 1.0))
    big_d = little_d / math.sqrt(vd)

    return big_d

# get variant site info
allele_freqs = {}

for line in pysam.VariantFile(vcf_file).fetch():
    chromo = line.contig

    alt_allele_freq = round(line.info['AF'][0], 2)
    ref_seq = line.ref
    alt_seq = line.alts[0]

    # if INDEL mode works out if polarised and if so if del or ins
    if mode == 'INDEL':
        variant_type = 'Unpolarised'
        try:
            ancestral_sequence = line.info['AA']
        except KeyError:
            ancestral_sequence = ''
    
        # determine type
        if len(alt_seq) == len(ancestral_sequence):
            if len(alt_seq) > len(ref_seq):
                variant_type = 'Del'
            elif len(alt_seq) < len(ref_seq):
                variant_type = 'Ins'
        elif len(ref_seq) == len(ancestral_sequence):
            if len(alt_seq) < len(ref_seq):
                variant_type = 'Del'
            elif len(alt_seq) > len(ref_seq):
                variant_type = 'Ins'

    # skips above for SNP mode
    else:
        variant_type = 'SNP'

    # determine region variant falls in
    if regions[0] != 'NONE':
        try:
            var_region = line.info['ANNO']
            if var_region.startswith('CDS'):
                var_region = 'CDS'
            if mode == 'SNP_comp':
                # catch S<->S and W<->W
                if var_region == 'intergenic':
                    if re.search(r'[CG]', alt_seq) and re.search(r'[CG]', ref_seq):
                        var_region = 'intergenic_ww_ss'
                    elif re.search(r'[AT]', alt_seq) and re.search(r'[AT]', ref_seq):
                        var_region = 'intergenic_ww_ss'
                    else:
                        continue

                # catch zerofold
                elif var_region == 'CDS':
                    try:
                        degen = line.info['DEGEN']
                        if degen == 0:
                            var_region = 'CDS_zerofold'
                        else:
                            continue
                    except KeyError:
                        continue

                # skip all others
                else:
                    continue
        except KeyError:
            continue
    else:
        var_region = 'NO_BINNING'

    # record allele freqs
    if chromo in allele_freqs.keys():
        if var_region in allele_freqs[chromo].keys():
            if variant_type in allele_freqs[chromo][var_region].keys():
                allele_freqs[chromo][var_region][variant_type].append(alt_allele_freq)
            else:
                allele_freqs[chromo][var_region][variant_type] = [alt_allele_freq]
        else:
            allele_freqs[chromo][var_region] = {variant_type: [alt_allele_freq]}
    else:
        allele_freqs[chromo] = {var_region: {variant_type: [alt_allele_freq]}}

    # record genome wide allele freqs
    if 'genome' in allele_freqs.keys():
        if var_region in allele_freqs['genome'].keys():
            if variant_type in allele_freqs['genome'][var_region].keys():
                allele_freqs['genome'][var_region][variant_type].append(alt_allele_freq)
            else:
                allele_freqs['genome'][var_region][variant_type] = [alt_allele_freq]
        else:
            allele_freqs['genome'][var_region] = {variant_type: [alt_allele_freq]}
    else:
        allele_freqs['genome'] = {var_region: {variant_type: [alt_allele_freq]}}

# set variables
tw_genome, pi_val_genome, tajD_genome = 0, 0, 0

# write header
if markdown is True:
    print '|region|bin|type|seg_sities|callable|theta_w|t_lwr|t_upr|pi|pi_lwr|pi_upr|tajD|tajD_lwr|tajD_upr' \
          '\n|:----|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|:---:|'
else:
    print 'region\tbin\ttype\tseg_sites\tcallable\ttheta_w\tt_lwr\tt_upr\tpi\tpi_lwr\tpi_upr\ttajD\ttajD_lwr\ttajD_upr'

# calc stats
type_list = []
genome_wide_callable_sites = {'NO_BINNING': [0, 0], 'intron': [0, 0], 'intergenic': [0, 0],
                              'CDS': [0, 0], 'CDS_zerofold': [0, 0], 'intergenic_ww_ss': [0, 0]}
for chromosome in allele_freqs.keys():
    if chromosome != 'genome':
        # get positons of CDS, introns and intergenic
        if database_path != 'NONE':
            chr_database = gffutils.FeatureDB(database_dict[chromosome], keep_order=True)
            chr_coord_dict = {}

            # get introns
            introns = chr_database.create_introns()
            merged_introns = chr_database.merge(introns, ignore_strand=True)
            chr_coord_dict['intron'] = [(x.start, x.end) for x in merged_introns]

            # get intergenic
            genes = chr_database.merge(chr_database.features_of_type('gene', order_by=('seqid', 'start', 'end')),
                                       ignore_strand=True)
            intergenic = chr_database.interfeatures(genes, new_featuretype='intergenic')
            merged_intergenic = chr_database.merge(intergenic, ignore_strand=True)
            chr_coord_dict['intergenic'] = [(x.start, x.end) for x in merged_intergenic]

            # get cds
            cds = chr_database.features_of_type('CDS', order_by=('seqid', 'start', 'end'))
            merged_cds = chr_database.merge(cds, ignore_strand=True)
            chr_coord_dict['CDS'] = [(x.start, x.end) for x in merged_cds]

        else:
            chr_coord_dict = {}

        for region in allele_freqs[chromosome].keys():

            callable_seq = call_fasta.fetch(chromosome)
            if region == 'NO_BINNING':
                callable_sites = callable_seq.count('2')

            else:
                callable_sites = sum([callable_seq[s[0]-1: s[1]].count('2')
                                      for s in chr_coord_dict[region.split('_')[0]]])

            genome_wide_callable_sites[region][0] += callable_sites

            # for indels gets no. unpolarised sites, sets to 0 for SNPs
            if mode == 'INDEL':  # todo fix correction
                unpoled = len(allele_freqs[chromosome][region]['Unpolarised'])
                type_list = ['Del', 'Ins', 'INDEL']
            else:
                unpoled = 0
                type_list = ['SNP']

            genome_wide_callable_sites[region][1] += unpoled

            # calculates statistics and prints to standard out
            for variant in type_list:
                if variant == 'INDEL':
                    freqs = sum([list(x) for x in allele_freqs[chromosome][region].values()], [])
                    unpoled = 0
                else:
                    try:
                        freqs = allele_freqs[chromosome][region][variant]
                    except KeyError:
                        if markdown is True:
                            print '|'+'|'.join([chromosome, region, variant, '0',
                                                str(callable_sites-unpoled), 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                                                'NA', 'NA', 'NA'])+'|'
                        else:
                            print '\t'.join([chromosome, region, variant, 'NA', 'NA', 'NA', 'NA', 'NA', 'NA',
                                             'NA', 'NA', 'NA'])
                        continue

                # bootstrapping

                if bootstrap == 0:
                    tw = theta_w(20, len(freqs))/(callable_sites-unpoled)
                    pi_val = pi(20, freqs)/(callable_sites-unpoled)
                    tajD = tajimas_d(20, freqs)

                    # set mean and cis for output if no bootsrapping
                    mean_tw = tw
                    ci_tw = [0, 0]
                    mean_pi = pi_val
                    ci_pi = [0, 0]
                    mean_tajD = tajD
                    ci_tajD = [0, 0]

                else:
                    bs_theta = []
                    bs_pi = []
                    bs_tajD = []

                    # resample with replacement - calc stats for each resample
                    for bs in range(0, bootstrap):
                        resampled_sfs = []
                        for i in range(0, len(freqs)):
                            random_no = random.randint(0, len(freqs)-1)
                            resampled_sfs.append(freqs[random_no])
                        tw = theta_w(20, len(resampled_sfs))/(callable_sites-unpoled)
                        bs_theta.append(tw)
                        pi_val = pi(20, resampled_sfs)/(callable_sites-unpoled)
                        bs_pi.append(pi_val)
                        tajD = tajimas_d(20, resampled_sfs)
                        bs_tajD.append(tajD)

                    # calc mean and cis
                    mean_tw = numpy.mean(bs_theta)
                    ci_tw = numpy.percentile(bs_theta, [2.5, 97.5])
                    mean_pi = numpy.mean(bs_pi)
                    ci_pi = numpy.percentile(bs_pi, [2.5, 97.5])
                    mean_tajD = numpy.mean(bs_tajD)
                    ci_tajD = numpy.percentile(bs_tajD, [2.5, 97.5])

                if markdown is True:
                    print '|'+'|'.join([chromosome, region, variant, str(len(freqs)), str(callable_sites-unpoled),
                                        str(round(mean_tw, 5)), str(round(ci_tw[0], 5)), str(round(ci_tw[1], 5)),
                                        str(round(mean_pi, 5)), str(round(ci_pi[0], 5)), str(round(ci_pi[1], 5)),
                                        str(round(mean_tajD, 3)), str(round(ci_tajD[0], 5)),
                                        str(round(ci_tajD[1], 5)),
                                        ])+'|'
                else:
                    print '\t'.join([chromosome, region, variant, str(len(freqs)), str(callable_sites-unpoled),
                                     str(mean_tw), str(ci_tw[0]), str(ci_tw[1]),
                                     str(mean_pi), str(ci_pi[0]), str(ci_pi[1]),
                                     str(mean_tajD), str(ci_tajD[0]), str(ci_tajD[1]),
                                     ])

# genome wide calculations
for region in allele_freqs['genome'].keys():

    # calculates statistics and prints to standard out
    for variant in type_list:
        if variant == 'INDEL':
            freqs = sum([list(x) for x in allele_freqs['genome'][region].values()], [])
            unpoled = 0
        else:
            freqs = allele_freqs['genome'][region][variant]
            unpoled = genome_wide_callable_sites[region][1]

        # bootstrapping
        if bootstrap == 0:
            tw = theta_w(20, len(freqs))/(genome_wide_callable_sites[region][0]-unpoled)
            pi_val = pi(20, freqs)/(genome_wide_callable_sites[region][0]-unpoled)
            tajD = tajimas_d(20, freqs)

            # set mean and cis for output if no bootsrapping
            mean_tw = tw
            ci_tw = [0, 0]
            mean_pi = pi_val
            ci_pi = [0, 0]
            mean_tajD = tajD
            ci_tajD = [0, 0]

        else:
            bs_theta = []
            bs_pi = []
            bs_tajD = []

            # resample with replacement - calc stats for each resample
            for bs in range(0, bootstrap):
                resampled_sfs = []
                for i in range(0, len(freqs)):
                    random_no = random.randint(0, len(freqs)-1)
                    resampled_sfs.append(freqs[random_no])
                tw = theta_w(20, len(resampled_sfs))/(genome_wide_callable_sites[region][0]-unpoled)
                bs_theta.append(tw)
                pi_val = pi(20, resampled_sfs)/(genome_wide_callable_sites[region][0]-unpoled)
                bs_pi.append(pi_val)
                tajD = tajimas_d(20, resampled_sfs)
                bs_tajD.append(tajD)

            # calc mean and cis
            mean_tw = numpy.mean(bs_theta)
            ci_tw = numpy.percentile(bs_theta, [2.5, 97.5])
            mean_pi = numpy.mean(bs_pi)
            ci_pi = numpy.percentile(bs_pi, [2.5, 97.5])
            mean_tajD = numpy.mean(bs_tajD)
            ci_tajD = numpy.percentile(bs_tajD, [2.5, 97.5])

        if markdown is True:
            print '|'+'|'.join(['genome', region, variant, str(len(freqs)),
                                str(genome_wide_callable_sites[region][0]-unpoled),
                                str(round(mean_tw, 5)), str(round(ci_tw[0], 5)), str(round(ci_tw[1], 5)),
                                str(round(mean_pi, 5)), str(round(ci_pi[0], 5)), str(round(ci_pi[1], 5)),
                                str(round(mean_tajD, 3)), str(round(ci_tajD[0], 5)),
                                str(round(ci_tajD[1], 5)),
                                ])+'|'
        else:
            print '\t'.join(['genome', region, variant, str(len(freqs)),
                             str(genome_wide_callable_sites[region][0]-unpoled),
                             str(mean_tw), str(ci_tw[0]), str(ci_tw[1]),
                             str(mean_pi), str(ci_pi[0]), str(ci_pi[1]),
                             str(mean_tajD), str(ci_tajD[0]), str(ci_tajD[1]),
                             ])
