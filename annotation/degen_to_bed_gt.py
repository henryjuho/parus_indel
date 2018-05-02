#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip
import os


def contig_dict(key_file):

    contigs = {}

    for line in open(key_file):

        refseq, chr_id = line.rstrip().split()

        contigs[refseq] = chr_id

    return contigs


def degeneracy(codon):

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

    ref_codon = codon[0].upper()
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


def cds_coord_to_codon_coords(fa_coord, direction):

    iter_coord = [x.split('..') for x in fa_coord.split(',')]

    # catch lone coords
    iter_coord_new = []
    for pos_pair in iter_coord:
        if len(pos_pair) == 1:
            iter_coord_new.append([pos_pair[0], pos_pair[0]])
        else:
            iter_coord_new.append(pos_pair)

    all_positions = []
    for block in iter_coord_new:
        block_poss = range(int(block[0]), int(block[1])+1)
        all_positions += block_poss
    if direction == 'complement':
        all_positions = all_positions[::-1]
    return all_positions


def pos_unique(pos, chr_degen_dict, other_degens):
    for z in other_degens:
        try:
            if pos in chr_degen_dict[z]:
                return False

        # some genes don't have all types of degen, particularly 3fold
        except KeyError:
            continue

    return True


def start_stop_ok(fa_seq):

    codon_list = [fa_seq[i:i+3] for i in range(0, len(fa_seq), 3)]
    stops = {'TAG', 'TGA', 'TAA'}
    # check if has start codon
    if codon_list[0] != 'ATG':
        return False

    # check if seq ends in stop codon
    elif codon_list[-1] not in stops:
        return False

    # check for premature stop
    else:
        for c in codon_list[1:-1]:
            if c in stops:
                return False

        return True


def main():
    # arguments
    parser = argparse.ArgumentParser(description='Script that outputs the position of '
                                                 'degenerate sites in coding regions to a bed file')
    parser.add_argument('-cds_fa', help='Fasta file with CDS sequences in', required=True)
    parser.add_argument('-degen', help='Degeneracy of sites to extract positions for', required=True,
                        action='append', choices=[0, 2, 3, 4], type=int)
    parser.add_argument('-contig_key', help=argparse.SUPPRESS,
                        default='{}/parus_indel/annotation/contigs_key.txt'.format(os.getenv('HOME')))
    args = parser.parse_args()

    # variables
    fa = args.cds_fa
    out_degens = set(args.degen)
    degen_data = {}
    contig_data = contig_dict(args.contig_key)

    # loop through fasta
    sequence, chromo, coords, skip, trans_name = '', '', [], False, ''
    for line in gzip.open(fa):

        # skip sequence lines following skipped header
        if skip is True:
            if not line.startswith('>'):
                continue
            else:
                skip = False

        # if a header line
        if line.startswith('>'):

            # process prev sequence
            if sequence != '' and len(sequence) % 3 == 0 and start_stop_ok(sequence):
                for i in range(0, len(sequence), 3):
                    codon = sequence[i:i+3]

                    # skip codons with Ns
                    if 'N' in codon:
                        continue

                    base_pos = coords[i:i+3]
                    for pos in range(0, 3):
                        degen = degeneracy([codon, pos])
                        site_pos = base_pos[pos]
                        if chromo not in degen_data.keys():
                            degen_data[chromo] = {}
                        if trans_name not in degen_data[chromo].keys():
                            degen_data[chromo][trans_name] = {}
                        if degen not in degen_data[chromo][trans_name].keys():
                            degen_data[chromo][trans_name][degen] = set()
                        degen_data[chromo][trans_name][degen] |= {site_pos}

            # skips records with partial cds or low quality
            if 'partial' in line or 'LOW QUALITY PROTEIN' in line:
                skip = True
                continue

            # reset holders and store details of next sequence
            sequence = ''
            ref_seq_chromo = '_'.join(line.split('[')[0].split('|')[1].split('_')[0:2])

            if ref_seq_chromo not in contig_data.keys():
                skip = True
                continue

            chromo = contig_data[ref_seq_chromo]
            trans_name = line.split('ID:')[1].split('[')[0].rstrip().rstrip(']')
            coord_data = line.split('[')[-1].rstrip().rstrip(']').replace('location=', '')

            # determine if complement or not
            if 'complement' in coord_data:
                method = 'complement'
            else:
                method = 'join'

            coord_string = coord_data.split('(')[-1].rstrip(')')
            coords = cds_coord_to_codon_coords(coord_string, method)

        # if a sequence line add to seq string
        else:
            sequence += line.rstrip()

    # process last seq in file
    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]

        # skip codons with Ns
        if 'N' in codon:
            continue

        base_pos = coords[i:i + 3]
        for pos in range(0, 3):
            degen = degeneracy([codon, pos])
            site_pos = base_pos[pos]
            if chromo not in degen_data.keys():
                degen_data[chromo] = {}
            if trans_name not in degen_data[chromo].keys():
                degen_data[chromo][trans_name] = {}
            if degen not in degen_data[chromo][trans_name].keys():
                degen_data[chromo][trans_name][degen] = set()
            degen_data[chromo][trans_name][degen] |= {site_pos}

    # output sites
    for contig in sorted(degen_data.keys()):
        for trans in degen_data[contig].keys():
            try:
                for degen_cat in [(degen_data[contig][trans][x], x) for x in out_degens]:
                    other_degen_cats = [0, 2, 3, 4]
                    other_degen_cats.remove(degen_cat[1])
                    positions = degen_cat[0]
                    for position in list(positions):
                        if pos_unique(position, degen_data[contig][trans], other_degen_cats):
                            print(contig.replace(' ', ''), str(int(position)-1), position, trans, sep='\t')
            # catch degens not found in gene, mainly 3 fold
            except KeyError:
                continue

if __name__ == '__main__':
    main()
