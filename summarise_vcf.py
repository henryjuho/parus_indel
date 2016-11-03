#!/usr/bin/env python

import argparse
from qsub import *
import sys
import math
import pysam

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='Vcf file to get summary stats for', required=True)
parser.add_argument('-call_fa', help='Fasta file of callable sites, coded as 0=N, 1=notcalled and 3=called',
                    required=True)
parser.add_argument('-mode', help='Variant mode, either SNP or INDEL', required=True)
parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
parser.add_argument('-md', help='If specified will output in markdown table format', default=False, action='store_true')
args = parser.parse_args()

# submission loop
if args.sub is True:
    command_line = [' '.join([y for y in sys.argv if y != '-sub'])]
    q_sub(command_line, 'set_output_here')
    sys.exit('Script submitted')

# variables
vcf_file = args.vcf
call_fasta = pysam.FastaFile(args.call_fa)
mode = args.mode
markdown = args.md


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

    # record allele freqs
    if chromo in allele_freqs.keys():
        if variant_type in allele_freqs[chromo].keys():
            allele_freqs[chromo][variant_type].append(alt_allele_freq)
        else:
            allele_freqs[chromo][variant_type] = [alt_allele_freq]
    else:
        allele_freqs[chromo] = {variant_type: [alt_allele_freq]}

# set variables
tw_genome, pi_val_genome, tajD_genome, genome_wide_callable_sites = 0, 0, 0, 0

# write header
if markdown is True:
    print '|region|type|theta_w|pi|tajD|\n|:----|:---:|:---:|:---:|:---:|'
else:
    print 'region\ttype\ttheta_w\tpi\ttajD'

# calc stats
for region in allele_freqs.keys():
    callable_sites = call_fasta.fetch(region).count('2')

    # for indels gets no. unpolarised sites, sets to 0 for SNPs
    if mode == 'INDEL':
        unpoled = len(allele_freqs[region]['Unpolarised'])
        type_list = ['Del', 'Ins', 'INDEL']
    else:
        unpoled = 0
        type_list = ['SNP']

    # calculates statistics and prints to standard out
    for variant in type_list:
        if variant == 'INDEL':
            freqs = sum([list(x) for x in allele_freqs[region].values()], [])
        else:
            try:
                freqs = allele_freqs[region][variant]
            except KeyError:
                if markdown is True:
                    print '|'+'|'.join([region, variant, 'NA', 'NA', 'NA'])+'|'
                else:
                    print '\t'.join([region, variant, 'NA', 'NA', 'NA'])
                continue

        tw = theta_w(20, len(freqs))/(callable_sites-unpoled)
        pi_val = pi(20, freqs)/(callable_sites-unpoled)
        tajD = tajimas_d(20, freqs)

        if markdown is True:
            print '|'+'|'.join([region, variant, str(tw), str(pi_val), str(tajD)])+'|'
        else:
            print '\t'.join([region, variant, str(tw), str(pi_val), str(tajD)])
