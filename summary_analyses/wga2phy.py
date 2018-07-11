#!/usr/bin/env python

from __future__ import print_function
import argparse
from qsub import q_sub
import sys
import gzip
import pysam
import os
sys.path.insert(0, os.getenv('HOME') + '/parus_indel/annotation')
from degen_to_bed_gt import cds_coord_to_codon_coords, start_stop_ok, contig_dict


def get_complement(dna_seq):

    dna_dict = {'A': 'T', 'a': 't', 'T': 'A', 't': 'a',
                'C': 'G', 'c': 'g', 'G': 'C', 'g': 'c',
                'N': 'N', 'n': 'n', '-': '-', '?': '?'}

    comp = ''.join([dna_dict[x] for x in dna_seq])

    return comp


def main():
    # arguments
    parser = argparse.ArgumentParser(description='Script that loops through a fasta file of CDS sequence and extracts '
                                                 'alignments for each sequence from a whole genome alignment bed file '
                                                 'and writes each alignment to a separate phylip file.')
    parser.add_argument('-wga_bed',
                        help='Whole genome alignment bed file to extract CDS alignments from',
                        required=True)
    parser.add_argument('-cds_fa',
                        help='Fasta file of CDS transcript sequences',
                        required=True)
    parser.add_argument('-out_dir',
                        help='Output directory to create and write phylip files to',
                        required=True)
    parser.add_argument('-sub',
                        help='If specified will submit script to cluster',
                        action='store_true',
                        default=False)
    parser.add_argument('-evolgen',
                        help='If specified will submit script to lab queue',
                        default=False,
                        action='store_true')
    parser.add_argument('-contig_key',
                        help=argparse.SUPPRESS,
                        default='{}/parus_indel/annotation/contigs_key.txt'.format(os.getenv('HOME')))
    args = parser.parse_args()

    # submission loop
    if args.sub is True:
        command_line = [' '.join([x for x in sys.argv if x != '-sub' and x != '-evolgen'])]
        q_sub(command_line, out=args.out_dir + 'cds_alignments_from_wga', evolgen=args.evolgen, t=96)
        sys.exit()

    # variables
    fa = args.cds_fa
    wga = pysam.TabixFile(args.wga_bed)
    out_dir = args.out_dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
    contig_data = contig_dict(args.contig_key)

    # loop through fasta
    sequence, chromo, coords, skip, trans_name, method = '', '', [], False, '', ''
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

                align_dict = {}

                # extract positions from alignment
                for codon in [coords[i: i + 3] for i in range(0, len(coords), 3)]:
                    codon_dict = {}
                    for pos in codon:
                        align_pos = list(wga.fetch(chromo, pos - 1, pos, parser=pysam.asTuple()))

                        # get align data in form [('dmel', 'T'), ('dsim', '?'), ('dyak', 'T')]

                        try:
                            spp_bases = zip(align_pos[0][4].split(','), align_pos[0][7].split(','))
                            for spp in spp_bases:
                                if spp[0] not in codon_dict.keys():
                                    codon_dict[spp[0]] = ''
                                if len(align_pos) == 0:
                                    codon_dict[spp[0]] += 'N'
                                else:
                                    codon_dict[spp[0]] += spp[1][0]
                        except IndexError:
                            continue

                    # if entire codon is not in alignment skip
                    if len(codon_dict.keys()) == 0:
                        continue

                    # if codon is less than 3bp ie if positions not in alignment for all spp skip
                    if len(codon_dict.values()[0]) != 3:
                        continue

                    # check codons and add to sequence string
                    codon_base_content = ''.join(codon_dict.values())
                    if '-' not in codon_base_content and 'N' not in codon_base_content:
                        for s in codon_dict.keys():
                            if s not in align_dict.keys():
                                align_dict[s] = ''
                            align_dict[s] += codon_dict[s]

                # complement (will have already been reversed) if necessary
                if method == 'complement':
                    align_dict = {x: get_complement(align_dict[x]) for x in align_dict.keys()}

                # skip if whole gene missing from alignment
                if len(align_dict.values()) == 0:
                    continue

                # output seq file
                seq_len = len(align_dict.values()[0])
                n_spp = len(align_dict.keys())
                file_name = '{}{}.{}.{}'.format(out_dir, trans_name, seq_len, 'phy')
                with open(file_name, 'w') as phy:
                    print('\t{}\t{}'.format(n_spp, seq_len), file=phy)
                    for spp_seq in align_dict.keys():
                        out_seq = align_dict[spp_seq]
                        line_wrap = 60
                        print(spp_seq, file=phy)
                        for i in range(0, len(out_seq), line_wrap):
                            print(out_seq[i: i + line_wrap], file=phy)

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
    if sequence != '' and len(sequence) % 3 == 0 and start_stop_ok(sequence):

        align_dict = {}

        # extract positions from alignment
        for codon in [coords[i: i + 3] for i in range(0, len(coords), 3)]:
            codon_dict = {}
            for pos in codon:
                align_pos = list(wga.fetch(chromo, pos - 1, pos, parser=pysam.asTuple()))

                # get align data in form [('dmel', 'T'), ('dsim', '?'), ('dyak', 'T')]

                try:
                    spp_bases = zip(align_pos[0][4].split(','), align_pos[0][7].split(','))
                    for spp in spp_bases:
                        if spp[0] not in codon_dict.keys():
                            codon_dict[spp[0]] = ''
                        if len(align_pos) == 0:
                            codon_dict[spp[0]] += 'N'
                        else:
                            codon_dict[spp[0]] += spp[1][0]
                except IndexError:
                    continue

            # if entire codon is not in alignment skip
            if len(codon_dict.keys()) == 0:
                continue

            # if codon is less than 3bp ie if positions not in alignment for all spp skip
            if len(codon_dict.values()[0]) != 3:
                continue

            # check codons and add to sequence string
            codon_base_content = ''.join(codon_dict.values())
            if '-' not in codon_base_content and 'N' not in codon_base_content:
                for s in codon_dict.keys():
                    if s not in align_dict.keys():
                        align_dict[s] = ''
                    align_dict[s] += codon_dict[s]

        # complement (will have already been reversed) if necessary
        if method == 'complement':
            align_dict = {x: get_complement(align_dict[x]) for x in align_dict.keys()}

        # skip if whole gene missing from alignment
        if len(align_dict.values()) != 0:

            # output seq file
            seq_len = len(align_dict.values()[0])
            n_spp = len(align_dict.keys())
            file_name = '{}{}.{}.{}'.format(out_dir, trans_name, seq_len, 'phy')
            with open(file_name, 'w') as phy:
                print('\t{}\t{}'.format(n_spp, seq_len), file=phy)
                for spp_seq in align_dict.keys():
                    out_seq = align_dict[spp_seq]
                    line_wrap = 60
                    print(spp_seq, file=phy)
                    for i in range(0, len(out_seq), line_wrap):
                        print(out_seq[i: i + line_wrap], file=phy)

if __name__ == '__main__':
    main()
