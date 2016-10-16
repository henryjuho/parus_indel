#!/usr/bin/env python

import argparse
from qsub import *
import sys
from pysam import VariantFile

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='Allsites vcf to apply filters to and get callable sites', required=True)
parser.add_argument('-bed', '--bed_repeats', help='BED file with repeat regions listed', required=True)
parser.add_argument('-DF', '--DepthFilter',
                    help='Defines abnormal depth eg) 2 means abnormal depth is twice and half the mean depth',
                    default=2.0, type=float)
parser.add_argument('-mean_depth', '--mean_depth', help='Mean coverage depth of samples', default=44.0)
parser.add_argument('-N', '--no_individuals', help='Number of individuals in allsites VCF', type=float, default=10.0)
parser.add_argument('-chr', help='Specifies chromosome to extract callable sites for, if ALL will run a job for each, '
                                 '-chr ALL can only be specified in conjunction with -sub', default='ALL')
parser.add_argument('-out', help='Output directory and prefix', required=True)
parser.add_argument('-evolgen', help='If specified will run on lab queue', action='store_true', default=False)
parser.add_argument('-sub', help='If specified will submit itself to cluster', action='store_true', default=False)
args = parser.parse_args()

# variables
all_sites = args.vcf
repeat_bed = args.bed_repeats
filter_factor = args.DepthFilter
all_data_mean_depth = float(args.mean_depth)
no_indiv = args.no_individuals
chromosome = args.chr
out = args.out
fasta_out = out + '.' + chromosome + '.fa'
evolgen = args.evolgen


# functions
def in_repeat(start, rep_list):
    for rep_region in rep_list:
        if rep_region[1] <= start <= rep_region[2]:
            return True
    else:
        return False

# submission loop
if args.sub is True:
    if chromosome == 'ALL':

        # gen chromo list and submit job for each
        grep_cmd = ('zcat ' + all_sites +
                    ' | head -n 4000 | grep ^##contig | cut -d "," -f 1 | cut -d "=" -f 3 | grep ^c')
        chromo_list = subprocess.Popen(grep_cmd, stdout=subprocess.PIPE, shell=True).communicate()[0].split('\n')[:-1]
        chromo_list = [x for x in chromo_list if x.startswith('chr')]
        output_fasta_list = []
        jid_list = []
        for chromo in chromo_list:
            output_fasta_list.append(out + '.' + chromo + '.fa')
            jid = chromo + '.callablesites.sh'
            jid_list.append(jid)
            command_line = ('./callable_sites_from_vcf.py '
                            '-vcf ' + all_sites + ' '
                            '-bed ' + repeat_bed + ' '
                            '-DF ' + str(filter_factor) + ' '
                            '-mean_depth ' + str(all_data_mean_depth) + ' '
                            '-N ' + str(no_indiv) + ' '
                            '-chr ' + chromo + ' '
                            '-out ' + out)
            q_sub([command_line], out + '.' + chromo, evolgen=evolgen, t=48)

        # cat job for final output
        cat_cmd = 'cat ' + ' '.join(output_fasta_list) + ' > ' + fasta_out
        q_sub([cat_cmd], out + 'cat', evolgen=evolgen, hold=jid_list)
        sys.exit()

    else:
        # submit script for chromosome
        command_line = ('./callable_sites_from_vcf.py '
                        '-vcf ' + all_sites + ' '
                        '-bed ' + repeat_bed + ' '
                        '-DF ' + str(filter_factor) + ' '
                        '-mean_depth ' + str(all_data_mean_depth) + ' '
                        '-N ' + str(no_indiv) + ' '
                        '-chr ' + chromosome + ' '
                        '-out ' + out)
        q_sub([command_line], out, evolgen=evolgen, t=48)
        sys.exit()

# catch -all specified without -sub
if args.chr == 'ALL' and args.sub is False:
    sys.exit('"-chr ALL" can only be run in conjunction with "-sub"')

# calculate depth cutoffs
lower_depth_limit = all_data_mean_depth / filter_factor
upper_depth_limit = all_data_mean_depth * filter_factor

# get bed regions per chromo
repeats = [(x.split()[0], x.split()[1], x.split()[2]) for x in open(repeat_bed) if x.split()[0] == chromosome]

# loop through allsites for chromosome
counter = 0
seq_count = 0
fasta_string = '>' + chromosome
with open(fasta_out, 'w') as out_fa:
    for line in VariantFile(all_sites).fetch(chromosome):

        # add line break every 60 bases
        if counter % 60 == 0:
            out_fa.write(fasta_string + '\n')
            fasta_string = ''
        counter += 1

        # check for ns
        if line.ref == 'N':
            fasta_string += '0'
            continue

        # depth filter
        try:
            cumulative_depth = line.info["DP"]
        except KeyError:
            fasta_string += '1'
            continue

        locus_mean_depth = cumulative_depth / no_indiv
        if lower_depth_limit <= locus_mean_depth <= upper_depth_limit:

            # repeat filter
            if in_repeat(line.pos, repeats) is False:
                fasta_string += '2'
                continue

            else:
                fasta_string += '1'
                continue

        else:
            fasta_string += '1'
            continue

    out_fa.write(fasta_string + '\n')

print counter
