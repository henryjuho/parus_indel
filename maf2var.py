#!/usr/bin/env python

import argparse
import subprocess

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-bed', help='Bed file with regions to extract', required=True)
parser.add_argument('-maf', help='Multiple alignment file to extract regions from', required=True)
args = parser.parse_args()

# variables
bed = args.bed
regions = [region.split() for region in open(bed)]
target_spp = regions[0][0].split('.')[0]
output_file = bed.replace('.bed', '.region.seq')
output = open(output_file, 'w')
maf_name = args.maf


# functions
def get_seq(maf_line):
    maf_line = maf_line.split()
    seqs = maf_line[1]
    return seqs


def get_start(maf_line):
    maf_line = maf_line.split()
    start = int(maf_line[2])
    return start


def get_end(maf_line):
    start = get_start(maf_line)
    maf_line = maf_line.split()
    end = start + int(maf_line[3])
    return end


def find_str_len(seqs, bp_length):
    base_count = 0
    matched_seq = ''
    for base in seqs:
        if base_count <= bp_length:
            matched_seq += base
            if base != '-':
                base_count += 1
    return len(matched_seq) - 1

# get species list
spp_grep = 'grep ^s ' + maf_name + ' | cut -d " " -f 2 | cut -d "." -f 1 | sort -u'
spp_list = sorted(subprocess.Popen(spp_grep, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1])

# write output header
header = 'ID\tstart\tend'
for spp in spp_list:
    header += '\t' + spp
header += '\n'
output.write(header)

# search for region in bed
s_data = []
a_data = ''

with open(maf_name) as maf:
    for line in maf:
        # extract from block
        if line.startswith('s'):
            s_data.append(line)

        # process block info
        else:
            if len(s_data) > 0:
                # identify position of target species in block
                position = -1
                target_pos = 'no_position'
                for spp in s_data:
                    position += 1
                    seq_spp = get_seq(spp).split('.')[0]
                    if seq_spp == target_spp:
                        target_pos = position

                # skip if spp not in block
                if target_pos == 'no_position':
                    s_data = []
                    a_data = ''
                    continue

                # get variant python coordinates
                target = str(s_data[target_pos])
                seq_ID = get_seq(target)
                start_pos = get_start(target)
                end_pos = get_end(target)
                for region in regions:
                    region_start = int(region[1])
                    region_end = int(region[2])

                    # identify regions which fall within alignment blocks
                    if region[0] == seq_ID and start_pos <= region_start < region_end < end_pos:

                        split_s_data = [s.rstrip('\n').split() for s in s_data]
                        dna_seq = split_s_data[target_pos][6]

                        # works out how far into maf block bed region starts
                        python_start = find_str_len(dna_seq, region_start - start_pos)

                        # get length of seq contained within the maf block
                        length = region_end - region_start
                        python_end = python_start + find_str_len(dna_seq[python_start:], length)
                        seq_data = {s[1].split('.')[0]:
                                    [s[6][python_start - 1], s[6][python_start: python_end], s[6][python_end]]
                                    for s in split_s_data}

                        output_line = '\t'.join(region)
                        for spp in spp_list:
                            try:
                                output_line += '\t' + ','.join(seq_data[spp])
                            except KeyError:
                                output_line += '\t.,.,.'
                        output_line += '\n'
                        output.write(output_line)

            # reset list before next block
            s_data = []
            a_data = ''

output.close()
