#!/usr/bin/env python

import argparse
import sys
from qsub import *

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='VCF file containing snps to polarise', required=True)
parser.add_argument('-align_data',
                    help='Text file in bed type format with alignment sequences corresponding to '
                         'SNP positions, generated using VARfromMAF.py',
                    required=True)
parser.add_argument('-target_spp', help='Species of samples in VCF file, corresponding to align data file header',
                    required=True)
parser.add_argument('-sub', help='If specified will submit script to cluster', action='store_true', default=False)
args = parser.parse_args()

# submission loop
if args.sub is True:
    command_line = [' '.join([x for x in sys.argv if x != '-sub'])]
    q_sub(command_line, args.vcf[:args.vcf.rfind('/') + 1] + 'snp_polarisation', evolgen=True, mem=10, rmem=10)
    sys.exit()

# variables
vcf_file = args.vcf
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
not_in_alignment = 0
no_ref_mismatch = 0
low_coverage = 0
ambiguous = 0

# loop through vcf file and annotate INDELs
previous_line = ''
for line in open(vcf_file):
    if line.startswith('#'):
        if line.startswith('##contig') and previous_line.startswith('##INFO'):
            new_info = '##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">\n'
            annotated_vcf.write(new_info)
        previous_line = line
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
            not_in_alignment += 1
            annotated_vcf.write(orig_line)
            continue

        # identify positions without full coverage
        neighbour_bp = [[x[0], x[2]] for x in alignment_info.values()]
        full_coverage = True
        for bp in neighbour_bp:
            if bp[0] == '.' or bp[1] == '.':
                full_coverage = False
                break

        # skips sites where ref allele differs from that in alignment
        if not ref == alignment_info[spp][1].rstrip('-').upper():
            no_ref_mismatch += 1
            annotated_vcf.write(orig_line)
            continue

        # Skip sites without full species coverage
        elif full_coverage is False:
            low_coverage += 1
            annotated_vcf.write(orig_line)
            continue

        else:
            # identify if ref or alt is ancestral
            target_seq = alignment_info[spp][1].upper()
            out_group_seqs = [alignment_info[out_spp][1].upper() for out_spp in alignment_info.keys() if out_spp != spp]
            ref_anc = True
            alt_anc = True

            # identify ref ancestral
            for sequence in out_group_seqs:
                if ref != sequence.rstrip('-'):
                    ref_anc = False
                    break

            # identify alt ancestral
            for sequence in out_group_seqs:
                if alt != sequence.rstrip('-'):
                    alt_anc = False
                    break

            # skip ambiguous sites
            if alt_anc is ref_anc:
                ambiguous += 1
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
                annotated_vcf.write(polarised_line)
            match += 1

print('\n'
      'Total no SNPs    : ' + str(counter) + '\n'
      'SNPs polarised   : ' + str(match) + '\n'
      'Ref mismatch     : ' + str(no_ref_mismatch) + '\n'
      'Low spp coverage : ' + str(low_coverage) + '\n'
      'Ambiguous        : ' + str(ambiguous) + '\n'
      'Not in alignment : ' + str(not_in_alignment) + '\n'
      'Total unpolarised: ' + str(counter - match) + '\n')

annotated_vcf.close()
