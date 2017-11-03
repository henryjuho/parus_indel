#!/usr/bin/env python

import argparse
import gffutils
import os
import gzip

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-gff', help='GFF annotation file', required=True)
parser.add_argument('-vcf', help='VCF file to annotate variants in', required=True)
parser.add_argument('-ref', help='Reference genome', required=True)
parser.add_argument('-chr', help='Chromosome to annotate, must be consistent with GFF and VCF', required=True)
parser.add_argument('-db_dir', help='Location of databases', required=True)
parser.add_argument('-out', help='Output directory', required=True)
args = parser.parse_args()

# variables
gff = args.gff
vcf = args.vcf
out_dir = args.out
db_dir = args.db_dir
ref = args.ref
chromo = args.chr
annotated_vcf = open(out_dir + vcf[vcf.rfind('/')+1:].replace('.vcf', '.degen.' + chromo + '.vcf'), 'w')
gff_db = db_dir + gff[gff.rfind('/')+1:].replace('.gff.gz', '.' + chromo + '.db')

# codon table
standard_codon_table = {
    "TTT": "F", "CTT": "L", "ATT": "I", "GTT": "V",
    "TTC": "F", "CTC": "L", "ATC": "I", "GTC": "V",
    "TTA": "L", "CTA": "L", "ATA": "I", "GTA": "V",
    "TTG": "L", "CTG": "L", "ATG": "M", "GTG": "V",

    "TCT": "S", "CCT": "P", "ACT": "T", "GCT": "A",
    "TCC": "S", "CCC": "P", "ACC": "T", "GCC": "A",
    "TCA": "S", "CCA": "P", "ACA": "T", "GCA": "A",
    "TCG": "S", "CCG": "P", "ACG": "T", "GCG": "A",

    "TAT": "Y", "CAT": "H", "AAT": "N", "GAT": "D",
    "TAC": "Y", "CAC": "H", "AAC": "N", "GAC": "D",
    "TAA": "*", "CAA": "Q", "AAA": "K", "GAA": "E",
    "TAG": "*", "CAG": "Q", "AAG": "K", "GAG": "E",

    "TGT": "C", "CGT": "R", "AGT": "S", "GGT": "G",
    "TGC": "C", "CGC": "R", "AGC": "S", "GGC": "G",
    "TGA": "*", "CGA": "R", "AGA": "R", "GGA": "G",
    "TGG": "W", "CGG": "R", "AGG": "R", "GGG": "G",
}


# functions
def snp_cds(snp_position, cds_coord_list):
    for cds_record in cds_coord_list:
        if cds_record[1] <= snp_position <= cds_record[2]:
            return cds_record

    return False


def reverse_comp(sequence):
    comp_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '-': '-', 'N': 'N'}
    sequence = sequence.upper()
    sequence_rev = ''
    for i in range(1, len(sequence)+1):
        sequence_rev += comp_dict[sequence[-i]]
    return sequence_rev


def snp_codon(snp_position, cds_record, reference):
    if cds_record[3] == '+':
        cds_seq = reference[cds_record[1]+cds_record[4]:cds_record[2]]
        snp_block_pos = snp_position - cds_record[1] - cds_record[4]

    else:
        cds_seq = reverse_comp(reference[cds_record[1]:cds_record[2]-cds_record[4]])
        snp_block_pos = cds_record[2]-1 - snp_position - cds_record[4]

    for i in range(0, len(cds_seq), 3):
        if i <= snp_block_pos < i+3:
            codon = cds_seq[i:i+3]
            snp_codon_pos = snp_block_pos - i
            return codon, snp_codon_pos


def degeneracy(codon):
    ref_codon = codon[0]
    snp_site = codon[1]
    ref_amino_acid = standard_codon_table[ref_codon]
    degen = 0
    for bp in ['A', 'T', 'C', 'G']:
        if snp_site == 0:
            alt_codon = bp + ref_codon[1:]
        elif snp_site == 1:
            alt_codon = ref_codon[0] + bp + ref_codon[2]
        else:
            alt_codon = ref_codon[0:2] + bp
        alt_amino_acid = standard_codon_table[alt_codon]
        if alt_amino_acid == ref_amino_acid:
            degen += 1
    if degen == 1:
        degen = 0
    return degen

# get chromosomal annotation data
gff_string = ''
for line in gzip.open(gff):
    gff_chromo = line.split('\t')[0]
    if chromo == gff_chromo:
        gff_string += line

# get chromosomal reference sequence
reference_sequence = ''
target = False
for line in open(ref):
    line = line.rstrip('\n')
    if line.startswith('>'):
        if line == '>' + chromo:
            target = True
        else:
            target = False
    else:
        if target is True:
            reference_sequence += line

# make / load gff database for chromo
if not os.path.isfile(gff_db):
    annotation_db = gffutils.create_db(gff_string, gff_db, merge_strategy='create_unique', from_string=True)
else:
    annotation_db = gffutils.FeatureDB(gff_db, keep_order=True)

# get cds
cds = annotation_db.features_of_type('CDS', order_by=('seqid', 'start', 'end'))
cds_intervals = list(set([(x.seqid, x.start-1, x.end, x.strand, int(x.frame)) for x in cds]))

# loop through vcf and identify if each variant is zerofold or not
previous_line = ''
all_variants = 0
cds_counter = 0
degen_counter = {0: 0, 2: 0, 3: 0, 4: 0}
for line in open(vcf):
    if line.startswith('#'):
        if line.startswith('##contig') and previous_line.startswith('##INFO'):
            new_info = '##INFO=<ID=DEGEN,Number=1,Type=Integer,Description="Annotation of SNP degeneracy">\n'
            annotated_vcf.write(new_info)
        previous_line = line
        annotated_vcf.write(line)
    elif line.split('\t')[0] == chromo:
        all_variants += 1
        split_line = line.split('\t')
        ref_allele = split_line[3].upper()
        alt_allele = split_line[4].upper()
        snp_pos = int(split_line[1])-1
        snp_cds_coords = snp_cds(snp_pos, cds_intervals)
        if snp_cds_coords is not False and split_line[7].find('ANNO=CDS_non_frameshift') != -1:
            cds_counter += 1
            snp_containing_codon = snp_codon(snp_pos, snp_cds_coords, reference_sequence)
            if snp_cds_coords[3] == '+':
                try:
                    assert snp_containing_codon[0][snp_containing_codon[1]] == ref_allele

                # catches instance where snp is in cds coords but due to gff 'frame' offset is not in exon
                except TypeError:
                    # write unnanotated line
                    annotated_vcf.write(line)
                    continue
            else:
                try:
                    assert snp_containing_codon[0][snp_containing_codon[1]] == reverse_comp(ref_allele)

                # catches instance where snp is in cds coords but due to gff 'frame' offset is not in exon
                except TypeError:
                    # write unnanotated line
                    annotated_vcf.write(line)
                    continue

            # catch out of frame cds
            if len(snp_containing_codon[0]) != 3:
                # write unnanotated line
                annotated_vcf.write(line)
                continue

            codon_degeneracy = degeneracy(snp_containing_codon)
            degen_counter[codon_degeneracy] += 1

            # write newly annotated line
            anno_line = ('\t'.join(split_line[0:7]) + '\t' + split_line[7] +
                         ';DEGEN=' + str(codon_degeneracy) + '\t' + '\t'.join(split_line[8:]))
            annotated_vcf.write(anno_line)

        else:
            # write unnanotated line
            annotated_vcf.write(line)

annotated_vcf.close()

print('\n'
      '|Category          | Number SNPs     |\n'
      '|:-----------------|:---------------:|\n'
      '|All               | ' + str(all_variants) + ' |\n'
      '|CDS               | ' + str(cds_counter) + ' |')

for z in degen_counter.keys():
    print'|' + str(z) + 'fold             | ' + str(degen_counter[z]) + ' |'
