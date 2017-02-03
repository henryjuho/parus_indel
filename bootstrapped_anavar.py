#!/usr/bin/env python

import argparse
from qsub import *
from hen_utils import *
import sys
import random
from collections import Counter
import time


# functions
def chunk_sfs(sfs, block_len):
    chunked_list = [sfs[i: i+block_len] for i in range(0, len(sfs), block_len)]
    return chunked_list


def block_bootstrap(chunked_list):
    resamp_sfs = []
    for i in range(0, len(chunked_list)):
        random_no = random.randint(0, len(chunked_list)-1)
        resamp_sfs += chunked_list[random_no]
    return resamp_sfs


def q_wait(cmd, no_job_limit):
    qstat_grep = "Qstat | grep -c ^' '"
    no_sub = subprocess.Popen(qstat_grep, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
    no_sub = int(no_sub)
    if no_sub >= no_job_limit:
        time.sleep(60)
        q_wait(cmd, no_job_limit)
    else:
        subprocess.call(cmd, shell=True)
    return

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-indel_vcf', '--indel_vcf',
                    help='The vcf file to generate bootstrapped INDEL SFS from',
                    required=True)
parser.add_argument('-snp_vcf', '--snp_vcf',
                    help='The vcf file to generate bootstrapped SNP SFS from',
                    required=True)
parser.add_argument('-neutral_type', '--neutral_type',
                    help='ID of neutral snps in snp sfs file, eg intergenic_ww_ss (default)',
                    default='intergenic_ww_ss')
parser.add_argument('-n', '--n', help='Sample size', required=True)
parser.add_argument('-r', '--region', help='Genomic region of INDELs', required=True,
                    choices=['CDS', 'CDS_non_frameshift', 'CDS_frameshift', 'intron', 'intergenic', 'no_bins'])
parser.add_argument('-r2', '--region2', help='Optional second region of INDELs to compare between spectra',
                    choices=['CDS', 'CDS_non_frameshift', 'CDS_frameshift', 'intron', 'intergenic', 'no_bins'],
                    default='none')
parser.add_argument('-lrt', '--likelihood_ratio_test', help='Paramter to fix for likelihood ratio test',
                    choices=['gamma_indel', 'gamma_ins', 'gamma_del', 'kappa', 'equal_gammas_ins', 'equal_gammas_del',
                             'equal_kappas'], required=True)
parser.add_argument('-out', '--output_prefix',
                    help='Output path and file prefix for control, data and results files',
                    required=True)
parser.add_argument('-n_boot', '--n_boot', help='Number of bootstrap replicates', required=True)
parser.add_argument('-block_size', '--block_size', help='Size of bootstrap block', default=100, type=int)
parser.add_argument('-sfs_path', '--sfs_path', help='Directory to hold bootstrapped sfs files', required=True)
parser.add_argument('-evolgen', '--evolgen', help='If specified will submit to lab queue', action='store_true',
                    default=False)
parser.add_argument('-node', '--node', help='Will run on specified nodes', default='0', type=str)
parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
args = parser.parse_args()


# input checks
if args.region2 != 'none' and not args.likelihood_ratio_test.startswith('equal'):
    sys.exit('if -region2 specified, -lrt must be one of equal_gammas_ins ,equal_gammas_del, equal_kappas')

# submission loop
if args.sub is True:
    command_line = [' '.join([x for x in sys.argv if x != '-sub'])]
    q_sub(command_line, args.output_prefix + '.controlscript', t=72)
    sys.exit()

# variables
indel_vcf = args.indel_vcf
snp_vcf = args.snp_vcf
region = args.region
region2 = args.region2
n = args.n
lrt = args.likelihood_ratio_test
out_prefix = args.output_prefix
node = args.node
n_boot = int(args.n_boot)
sfs_path = args.sfs_path
if not os.path.isdir(sfs_path):
    os.makedirs(sfs_path)
if not os.path.isdir(file_location(out_prefix)):
    os.makedirs(file_location(out_prefix))
evolgen = args.evolgen
block_size = args.block_size
pos_biallelic_freqs = [i/float(n) for i in range(1, int(n))]

# get sfs in
sfs_dict = {'snp': {}, 'del': {}, 'ins': {}}
snp_sfs_cmd = ('vcf2raw_sfs.py -vcf ' + snp_vcf + ' -region intergenic '
                                                  '-mode snp -mute_type SS -mute_type WW -auto_only')
snp_sfs = subprocess.Popen(snp_sfs_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
sfs_dict['snp']['intergenic_ww_ss'] = chunk_sfs(snp_sfs, block_size)

del_sfs_cmd = ('vcf2raw_sfs.py -vcf ' + indel_vcf + ' -region ' + region + ' -mode del -auto_only')
del_sfs = subprocess.Popen(del_sfs_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
sfs_dict['del'][region] = chunk_sfs(del_sfs, block_size)

ins_sfs_cmd = ('vcf2raw_sfs.py -vcf ' + indel_vcf + ' -region ' + region + ' -mode ins -auto_only')
ins_sfs = subprocess.Popen(ins_sfs_cmd, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
sfs_dict['ins'][region] = chunk_sfs(ins_sfs, block_size)

# second region
if region2 != 'none':
    del_sfs_cmd2 = ('vcf2raw_sfs.py -vcf ' + indel_vcf + ' -region ' + region2 + ' -mode del -auto_only')
    del_sfs2 = subprocess.Popen(del_sfs_cmd2, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    sfs_dict['del'][region2] = chunk_sfs(del_sfs2, block_size)

    ins_sfs_cmd2 = ('vcf2raw_sfs.py -vcf ' + indel_vcf + ' -region ' + region2 + ' -mode ins -auto_only')
    ins_sfs2 = subprocess.Popen(ins_sfs_cmd2, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    sfs_dict['ins'][region2] = chunk_sfs(ins_sfs2, block_size)

# write bootstrapped sfs files Freq\tProportion\tBin\tNorm
for replicate in range(0, n_boot+1):
    for var_type in sfs_dict.keys():
        sfs_output_file = sfs_path + var_type + '_rep' + str(replicate) + '.sfs.txt'
        with open(sfs_output_file, 'w') as sfs_out:
            sfs_out.write('Freq\tProportion\tBin\tNorm\n')
            for regional_sfs in sfs_dict[var_type].keys():
                if replicate == 0:
                    resampled_sfs = sum(sfs_dict[var_type][regional_sfs], [])
                else:
                    resampled_sfs = block_bootstrap(sfs_dict[var_type][regional_sfs])
                counted_sorted_sfs = sorted(Counter(resampled_sfs).most_common(), key=lambda z: z[0])
                sfs_freq_dict = {x[0]: x[1] for x in counted_sorted_sfs}

                for frequency in pos_biallelic_freqs:
                    try:
                        no_indels = sfs_freq_dict[str(frequency)]
                    except KeyError:
                        no_indels = 0  # account for freqs with 0 variants

                    sfs_out.write('\t'.join([str(frequency), str(no_indels), regional_sfs, 'NA']) + '\n')

    # write and submit anavar job
    anavar_cmd = ('anavar.py '
                  '-i_sfs ' + sfs_path + 'ins_rep' + str(replicate) + '.sfs.txt '
                  '-d_sfs ' + sfs_path + 'del_rep' + str(replicate) + '.sfs.txt '
                  '-s_sfs ' + sfs_path + 'snp_rep' + str(replicate) + '.sfs.txt '
                  '-n ' + n + ' '
                  '-r ' + region + ' '
                  '-lrt ' + lrt + ' '
                  '-out ' + out_prefix + '.' + region + '.' + lrt + '.rep' + str(replicate) + ' '
                  '-node ' + node + ' ')
    if region2 != 'none':
        anavar_cmd += '-r2 ' + region2 + ' '
    if evolgen is True:
        anavar_cmd += '-evolgen'

    # check how many jobs are submitted
    q_wait(anavar_cmd, 1000)

# write submission script for file merging
merge_cmd = 'merge_maxL.py -maxL_dir ' + file_location(out_prefix) + ' > ' + out_prefix + '.merged.maxL.txt'
q_write([merge_cmd], out=out_prefix+'.merge')
