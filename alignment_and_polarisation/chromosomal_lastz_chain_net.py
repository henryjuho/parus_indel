#!/usr/bin/env python

import argparse
import os
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-target_list',
                    help='List file (tab delim) containing fasta location and fasta nickname for target',
                    required=True)
parser.add_argument('-query_list',
                    help='List file (tab delim) containing fasta location and fasta nickname for query',
                    required=True)
parser.add_argument('-out',
                    help='Output directory',
                    required=True)
parser.add_argument('-evolgen',
                    help='Will run on evolgen queue',
                    action='store_true',
                    default=False)
args = parser.parse_args()

# variables
out_dir = args.out
if not os.path.isdir(out_dir):
    os.mkdir(out_dir)
evolgen = args.evolgen
target_genomes = [fa.rstrip('\n').split() for fa in open(args.target_list)]
query_genomes = [fa.rstrip('\n').split() for fa in open(args.query_list)]
alignment_list = []

# make alignment table
for target in target_genomes:
    for query in query_genomes:
        target_id = target[1].split('.')[1]
        query_id = query[1].split('.')[1]
        if target_id == query_id or target_id + 'A' == query_id or target_id + 'B' == query_id or \
           target_id == query_id + 'A' or target_id == query_id + 'B':

            alignment_list.append(target + query)

# iterate through alignments
for alignment in alignment_list:

    # identify target data and create 2bit files
    target = alignment[0]
    target_name = alignment[1]
    target_size_file = out_dir + target_name + '.sizes'
    if not os.path.isfile(target_size_file):
        subprocess.call('faSize -detailed ' + target + ' > ' + target_size_file, shell=True)
    target_2bit = out_dir + target_name + '.2bit'
    if not os.path.isfile(target_2bit):
        subprocess.call('faToTwoBit ' + target + ' ' + target_2bit, shell=True)

    # identify query data and create 2bit files
    query = alignment[2]
    query_name = alignment[3]
    query_size_file = out_dir + query_name + '.sizes'
    if not os.path.isfile(query_size_file):
        subprocess.call('faSize -detailed ' + query + ' > ' + query_size_file, shell=True)
    query_2bit = out_dir + query_name + '.2bit'
    if not os.path.isfile(query_2bit):
        subprocess.call('faToTwoBit ' + query + ' ' + query_2bit, shell=True)

    # determine output prefix
    out_prefix = out_dir + target_name + '.' + query_name

    # lastz commandline
    lastz_output = out_prefix + '.axt'
    lastz_cmd = ('lastz_latest ' +
                 target_2bit + '[nameparse=darkspace] ' +
                 query_2bit + '[nameparse=darkspace] '
                 '--format=axt '
                 '--step=19 '
                 '--hspthresh=2200 '
                 '--inner=2000 '
                 '--ydrop=3400 '
                 '--gappedthresh=10000 '
                 '--scores=/fastdata/bop15hjb/bird_alignments/UCSC_pipeline/Scores/HoxD55 '
                 '--chain '
                 '> ' + lastz_output)

    # chaining commandlines
    chain_out = out_prefix + '.chain'
    axtChain_cmd = ('axtChain '
                    '-linearGap=loose ' +
                    lastz_output + ' ' +
                    target_2bit + ' ' +
                    query_2bit + ' ' +
                    chain_out)

    prenet_out = out_prefix + '.prenet.chain'
    chainPreNet_cmd = ('chainPreNet ' +
                       chain_out + ' ' +
                       target_size_file + ' ' +
                       query_size_file + ' ' +
                       prenet_out)

    # netting commandline
    target_net_out = out_prefix + '.target.net'
    query_net_out = out_prefix + '.query.net'
    chain_net_cmd = ('chainNet ' +
                     prenet_out + ' ' +
                     target_size_file + ' ' +
                     query_size_file + ' ' +
                     target_net_out + ' ' +
                     query_net_out)

    syn_target_net_out = out_prefix + '.target.syn.net'
    syn_query_net_out = out_prefix + '.query.syn.net'
    target_netsyntenic_cmd = ('netSyntenic ' + target_net_out + ' ' + syn_target_net_out)
    query_netsyntenic_cmd = ('netSyntenic ' + query_net_out + ' ' + syn_query_net_out)

    # converting to maf format for use with multiz
    post_net_axt = out_prefix + '.postnet.axt'
    netToAxt_cmd = ('netToAxt ' +
                    syn_target_net_out + ' ' +
                    prenet_out + ' ' +
                    target_2bit + ' ' +
                    query_2bit + ' ' +
                    post_net_axt)

    sorting_out = out_prefix + '.postnet.sorted.axt'
    axtSort_cmd = ('axtSort ' + post_net_axt + ' ' + sorting_out)

    output_maf = out_prefix + '.maf'
    axtToMaf_cmd = ('axtToMaf ' +
                    sorting_out + ' ' +
                    target_size_file + ' ' +
                    query_size_file + ' ' +
                    output_maf)

    # submit to SGE
    if evolgen is True:
        q_sub([lastz_cmd, axtChain_cmd, chainPreNet_cmd, chain_net_cmd, target_netsyntenic_cmd,
              query_netsyntenic_cmd, netToAxt_cmd, axtSort_cmd, axtToMaf_cmd],
              out_prefix,
              t=72,
              rmem=10, mem=10,
              evolgen=True)

    else:
        q_sub([lastz_cmd, axtChain_cmd, chainPreNet_cmd, chain_net_cmd, target_netsyntenic_cmd,
              query_netsyntenic_cmd, netToAxt_cmd, axtSort_cmd, axtToMaf_cmd],
              out_prefix,
              t=72,
              rmem=10, mem=10)
