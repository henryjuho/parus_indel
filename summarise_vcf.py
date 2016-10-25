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
parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
args = parser.parse_args()

# submission loop
if args.sub is True:
    command_line = [' '.join([y for y in sys.argv if y != '-sub'])]
    q_sub(command_line, 'set_output_here')
    sys.exit('Script submitted')

# variables
vcf_file = args.vcf
call_fasta = pysam.FastaFile(args.call_fa)


# functions
def theta_w(n, segsites):
    # theta = S/a
    alpha = sum(1.0/x for x in range(1, n))
    theta = float(segsites) / alpha
    return theta


def pi(n, allele_frq_list):
    no_seqs_dif = [(1.0 - raf**2 - (1.0-raf)**2) * (n/(n-1.0)) for raf in allele_frq_list]
    seq_pi = sum(no_seqs_dif)
    return seq_pi


def tajimas_d(n, allele_frq_list):
    segsites = float(len(allele_frq_list))
    little_d = pi(n, allele_frq_list) - theta_w(n, segsites)

    a1 = sum(1.0/x for x in range(1, n))
    a2 = sum(1.0/x**2 for x in range(1, n))

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
no_segsites = {'genome': 0}
allele_freqs = {'genome': []}

for line in pysam.VariantFile(vcf_file).fetch():
    chromo = line.contig
    alt_allele_freq = round(line.info['AF'][0], 2)
    no_segsites['genome'] += 1
    allele_freqs['genome'].append(alt_allele_freq)
    try:
        no_segsites[chromo] += 1
        allele_freqs[chromo].append(alt_allele_freq)
    except KeyError:
        no_segsites[chromo] = 1
        allele_freqs[chromo] = [alt_allele_freq]

# Calc stats
tw_genome, pi_val_genome, tajD_genome, genome_wide_callable_sites = 0, 0, 0, 0

print '|region|theta_w|pi|tajD|\n|:----|:---:|:---:|:---:|'
for region in allele_freqs.keys():
    if region == 'genome':
        freqs = allele_freqs[region]
        tw_genome = theta_w(20, len(freqs))
        pi_val_genome = pi(20, freqs)
        tajD_genome = tajimas_d(20, freqs)
    else:
        callable_sites = call_fasta.fetch(region).count('2')
        genome_wide_callable_sites += callable_sites
        freqs = allele_freqs[region]
        tw = theta_w(20, len(freqs))/callable_sites
        pi_val = pi(20, freqs)/callable_sites
        tajD = tajimas_d(20, freqs)
        print '|'+'|'.join([region, str(tw), str(pi_val), str(tajD)])+'|'

print '|'+'|'.join(['genome', str(tw_genome/genome_wide_callable_sites),
                   str(pi_val_genome/genome_wide_callable_sites), str(tajD_genome)])+'|'
