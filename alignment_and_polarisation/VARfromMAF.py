#!/usr/bin/env python

import argparse
import os

from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF file to trawl maf for its variants', required=True)
parser.add_argument('-maf', help='MAF file with alignment to find variants in', required=True)
parser.add_argument('-target_spp', help='Name of target species as it appears in the maf file', required=True)
parser.add_argument('-out', help='Output directory', required=True)
parser.add_argument('-no_jobs', help='Number of jobs to split maf variant extraction over', required=True)
args = parser.parse_args()

# variables
vcf_name = args.vcf
maf = args.maf
target_spp = args.target_spp
output_dir = args.out
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)
bed_dir = output_dir + 'bed_files/'
if not os.path.isdir(bed_dir):
    os.mkdir(bed_dir)
bed_prefix = bed_dir + target_spp + '_bin_'
no_jobs = args.no_jobs

# calc number of variants per bed
variants_grep = 'grep -c -v ^# ' + vcf_name
no_vcf_variants = float(subprocess.Popen(variants_grep, shell=True, stdout=subprocess.PIPE).communicate()[0])
variants_per_bed = int(no_vcf_variants / float(no_jobs)) + 1

# vcf chr order
chr_order_grep = 'grep -v ^# ' + vcf_name + ' | cut -f 1 | uniq'
chr_list = subprocess.Popen(chr_order_grep, shell=True, stdout=subprocess.PIPE).communicate()[0].split()

# write beds
variant_tracker = 0
output_tracker = 1
bed_name = bed_prefix + str(output_tracker) + '.bed'
output_list = [bed_name.replace('.bed', '.region.seq')]
output_bed = open(bed_name, 'w')
bed_list = [bed_name]

with open(vcf_name) as vcf:
    for line in vcf:
        if not line.startswith('#'):
            # extract info from VCF line
            line = line.split()
            variant_tracker += 1
            seq_ID = target_spp + '.' + line[0]
            start = int(line[1]) - 1
            ref_len = len(line[3])
            end = str(start + ref_len)

            # construct bed line
            bed_line = '\t'.join([seq_ID, str(start), end]) + '\n'

            # write to beds
            if variant_tracker <= variants_per_bed:
                output_bed.write(bed_line)
            else:
                output_bed.close()
                output_tracker += 1
                variant_tracker = 1
                bed_name = bed_prefix + str(output_tracker) + '.bed'
                output_list.append(bed_name.replace('.bed', '.region.seq'))
                bed_list.append(bed_name)
                output_bed = open(bed_name, 'w')
                output_bed.write(bed_line)

# write and submit array job for extract_maf_region.py
extract_maf_region_cmd = ('~/parus_indel/alignment_and_polarisation/maf2var.py '
                          '-bed ' + bed_prefix + '$SGE_TASK_ID.bed '
                          '-maf ' + maf)

q_sub([extract_maf_region_cmd],
      bed_dir+'var_from_maf',
      t=72,
      jid='var_from_maf.sh',
      array=[1, int(no_jobs)])

# write hold job to concatonate all output files
concat_cmd = ('~/parus_indel/alignment_and_polarisation/concat_seq_files.py '
              '-out ' + output_dir + 'all_variants.alignment_states.txt ')

for output in output_list:
    concat_cmd += '-seq ' + output + ' '

q_sub([concat_cmd],
      output_dir + 'merge_alignment_variants',
      hold=['var_from_maf.sh'])
