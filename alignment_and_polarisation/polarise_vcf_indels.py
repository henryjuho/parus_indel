#!/usr/bin/env python

from __future__ import print_function
import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF file containing insertion and deletion variants to polarise', required=True)
parser.add_argument('-align_data',
                    help='Text file in bed type format with alignment sequences corresponding to '
                         'INDEL positions, generated using INDELsfromMAF.py',
                    required=True)
parser.add_argument('-target_spp', help='Species of samples in VCF file, corresponding to align data file header',
                    required=True)
parser.add_argument('-no_vcf', default=False, action='store_true', help=argparse.SUPPRESS)
args = parser.parse_args()

# variables
vcf_file = args.vcf
if not args.no_vcf:
    annotated_vcf = open(vcf_file.replace('.vcf', '.polarised.vcf'), 'w')
align_file = args.align_data
spp = args.target_spp

# build alignment data dictionary
spp_list = []
align_data = {}
for line in open(align_file):
    line = line.rstrip('\n').split('\t')
    if line[0] == 'ID':
        spp_list = line[3:]
    else:
        chromo = line[0].split('.')[1]
        vcf_position = str(int(line[1]) + 1)
        seq_key = chromo + '_' + vcf_position
        sequences = line[3:]
        align_data[seq_key] = {spp_list[i]: sequences[i].split(',') for i in range(0, len(spp_list))}

# counters
counter = 0
match = 0
no_hotspot = 0
low_coverage = 0
ambiguous = 0
not_aligned = 0

# loop through vcf file and annotate INDELs
previous_line = ''
for line in open(vcf_file):
    if line.startswith('#'):
        if line.startswith('##contig') and previous_line.startswith('##INFO'):
            new_info = '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\n'
            if not args.no_vcf:
                annotated_vcf.write(new_info)
        previous_line = line
        if not args.no_vcf:
            annotated_vcf.write(line)
    else:
        counter += 1
        orig_line = line
        line = line.split('\t')
        chrom = line[0]
        pos = line[1]
        ref = line[3]
        alt = line[4]
        info = line[7]
        try:
            alignment_info = align_data[chrom + '_' + pos]
        except KeyError:
            not_aligned += 1
            if not args.no_vcf:
                annotated_vcf.write(orig_line)
            continue

        # identify flanking INDELs and positions without full coverage
        neighbour_bp = [[x[0], x[2]] for x in alignment_info.values()]
        full_coverage = True
        neighbouring_deletions = False
        for bp in neighbour_bp:
            if bp[0] == '-' or bp[1] == '-':
                neighbouring_deletions = True
                break
            if bp[0] == '.' or bp[1] == '.':
                full_coverage = False
                break

        # identify indel hotspots
        indel_sequences = [y[1].rstrip('-') for y in alignment_info.values()]
        hotspot = False
        for indel in indel_sequences:
            if '-' in indel:
                hotspot = True
                break

        # skips sites where ref allele differs from that in alignment, ie insertion within INDEL
        if not ref == alignment_info[spp][1].rstrip('-').upper():
            no_hotspot += 1
            if not args.no_vcf:
                annotated_vcf.write(orig_line)
            continue

        # skip sites with multiple INDELs in outgroups
        elif hotspot is True:
            no_hotspot += 1
            if not args.no_vcf:
                annotated_vcf.write(orig_line)
            continue

        # Skip sites without full species coverage
        elif full_coverage is False:
            low_coverage += 1
            if not args.no_vcf:
                annotated_vcf.write(orig_line)
            continue

        # skips sites that are flanked by INDELs
        elif neighbouring_deletions is True:
            no_hotspot += 1
            if not args.no_vcf:
                annotated_vcf.write(orig_line)
            continue

        else:
            # identify if ref or alt is ancestral
            target_seq = alignment_info[spp][1]
            out_group_seqs = [alignment_info[out_spp][1] for out_spp in alignment_info.keys() if out_spp != spp]
            ref_anc = True
            alt_anc = True

            # identify ref ancestral
            for sequence in out_group_seqs:
                if len(ref) != len(sequence.rstrip('-')):
                    ref_anc = False
                    break

            # identify alt ancestral
            for sequence in out_group_seqs:
                if len(alt) != len(sequence.rstrip('-')):
                    alt_anc = False
                    break

            # skip ambiguous sites
            if alt_anc is ref_anc:
                ambiguous += 1
                if not args.no_vcf:
                    annotated_vcf.write(orig_line)
                continue

            # set AA annotation
            aa = 'NONE'
            if ref_anc is True:
                aa = ';AA=' + ref
            elif alt_anc is True:
                aa = ';AA=' + alt

            # write annotation
            if not aa == 'NONE':
                polarised_line = '\t'.join(line[0:7]) + '\t' + info + aa + '\t' + '\t'.join(line[8:])
                if not args.no_vcf:
                    annotated_vcf.write(polarised_line)
            match += 1

print('category\tcount\n'
      'total\t' + str(counter) + '\n'
      'polarised\t' + str(match) + '\n'
      'hotspot\t' + str(no_hotspot) + '\n'
      'low_coverage\t' + str(low_coverage) + '\n'
      'not_aligned\t' + str(not_aligned) + '\n'
      'ambiguous\t' + str(ambiguous) + '\n'
      'unpolarised\t' + str(counter - match) + '\n')

if not args.no_vcf:
    annotated_vcf.close()
