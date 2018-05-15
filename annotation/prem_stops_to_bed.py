#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip
import subprocess
import pysam
import os
from degen_to_bed_gt import start_stop_ok, cds_coord_to_codon_coords, contig_dict


def reverse_strand(codon_coords):

    """
    identifies if codon coordinates go backwards
    :param codon_coords: list
    :return: bool
    """

    if codon_coords[0] > codon_coords[1] > codon_coords[2]:
        return True

    elif codon_coords[0] < codon_coords[1] < codon_coords[2]:
        return False

    else:
        raise ValueError


def multi_nucleotide_poly(snp_position, snp_set, codon_pos):

    """
    identifies if adjacent bases in the codon are polymorphic
    :param snp_position: int
    :param snp_set: set
    :param codon_pos: int
    :return: bool
    """

    if codon_pos == 0 and snp_position + 1 in snp_set:
        return True

    elif codon_pos == 1:
        if snp_position + 1 in snp_set:
            return True
        if snp_position - 1 in snp_set:
            return True

    elif codon_pos == 2 and snp_position - 1 in snp_set:
        return True

    else:
        return False


def snp_makes_stop(pysam_snp, codon_pos, codon, codon_coords):

    """
    takes a biallelic snp and checks to see if it
    makes a premature stop codon
    :param pysam_snp: pysam.VariantRecord
    :param codon_pos: int
    :param codon: str
    :param codon_coords: list
    :return: bool
    """

    ref = pysam_snp.ref.upper()
    alt = pysam_snp.alts[0].upper()

    # reverse comp if neccessary
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    if reverse_strand(codon_coords):
        ref = comp[ref]
        alt = comp[alt]

    stop = prem_stop(codon, codon_pos, bases=(ref, alt))

    return stop


def prem_stop(codon, pos, bases=('A', 'T', 'G', 'C')):

    """
    goes through a list of bases and sees if any make a stop codon
    :param codon: str
    :param pos: int
    :param bases: set
    :return: bool
    """

    stops = {'TAG', 'TGA', 'TAA'}

    for base in bases:

        # construct possible codons
        if pos == 0:
            alt_codon = base + codon[1:]
        elif pos == 1:
            alt_codon = codon[0] + base + codon[2]
        elif pos == 2:
            alt_codon = codon[0:2] + base
        else:
            raise IndexError

        if alt_codon in stops:
            return True

    return False


def cds_snp_coords(vcf):

    """
    greps CDS snp coords from vcf and makes a dict of sets
    :param vcf: str
    :return: dict
    """

    grep = 'zgrep ANNO=CDS {} | cut -f 1,2'.format(vcf)

    cds_snps = subprocess.Popen(grep, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[:-1]
    cds_snps_dict = {}

    for snp in cds_snps:
        chromo, pos = snp.split()[0], int(snp.split()[1])
        if chromo not in cds_snps_dict.keys():
            cds_snps_dict[chromo] = set()

        cds_snps_dict[chromo] |= {pos}

    return cds_snps_dict


def process_transcript(sequence, coords, trans_name, nonsense_data, snps, chromo, vcf_records):

    """
    gets positions of possible nonsense mutation sites and nonsense snps for a transcript
    :param sequence: str
    :param coords: list
    :param trans_name: str
    :param nonsense_data: dict
    :param snps: set
    :param chromo: str
    :param vcf_records: pysam.VariantFile
    :return: None
    """

    if sequence != '' and len(sequence) % 3 == 0 and start_stop_ok(sequence):
        for i in range(0, len(sequence) - 3, 3):

            length = len(sequence)
            codon = sequence[i:i + 3]

            # skip codons with Ns
            if 'N' in codon:
                continue

            base_pos = coords[i:i + 3]

            for pos in range(0, 3):
                site_pos = base_pos[pos]
                if trans_name not in nonsense_data.keys():
                    nonsense_data[trans_name] = {'length': length, 'chromo': chromo, 'call': set(), 'snps': set()}

                # is a premature stop possible?
                if prem_stop(codon, pos):

                    nonsense_data[trans_name]['call'] |= {site_pos}

                    # if snp at site
                    if site_pos in snps[chromo]:

                        # check not an MNP
                        if multi_nucleotide_poly(site_pos, snps[chromo], pos):
                            continue

                        # does it introduce a stop?
                        snp_data = list(vcf_records.fetch(chromo, site_pos - 1, site_pos))
                        snp_record = snp_data[0]

                        # skip snps that don't make stop
                        if not snp_makes_stop(snp_record, pos, codon, base_pos):
                            continue

                        nonsense_data[trans_name]['snps'] |= {site_pos}

                    # move on if no snp at position
                    else:
                        continue

                # if site cannot give a premature stop
                else:
                    continue


def main():
    # arguments
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-cds_fa', help='Fasta file with CDS sequences in', required=True)
    parser.add_argument('-chr', help='chromosome', required=True)
    parser.add_argument('-vcf', help='SNP vcf path', required=True)
    parser.add_argument('-out', help='output file stem', required=True)
    parser.add_argument('-contig_key', help=argparse.SUPPRESS,
                        default='{}/parus_indel/annotation/contigs_key.txt'.format(os.getenv('HOME')))
    args = parser.parse_args()

    # variables
    fa = args.cds_fa
    nonsense_data = {}
    snps = cds_snp_coords(args.vcf)
    vcf_records = pysam.VariantFile(args.vcf)
    contig_data = contig_dict(args.contig_key)
    possible_nonsense = open('{}_{}_nonsense_potential.bed'.format(args.out, args.chr), 'w')
    actual_nonsense = open('{}_{}_nonsense_actual.bed'.format(args.out, args.chr), 'w')

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
            process_transcript(sequence, coords, trans_name, nonsense_data, snps, chromo, vcf_records)

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

            # skip non target chromos
            if chromo != args.chr:
                skip = True
                continue

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

    # process final sequence
    process_transcript(sequence, coords, trans_name, nonsense_data, snps, chromo, vcf_records)

    # output sites
    # with open(args.out, 'w') as out:
    for trans in sorted(nonsense_data.keys()):

        # trans_len = nonsense_data[trans]['length']
        nonse_chromo = nonsense_data[trans]['chromo']
        nonse_pos = nonsense_data[trans]['call']
        nonse_snp_pos = nonsense_data[trans]['snps']

        for positions in [[nonse_pos, possible_nonsense], [nonse_snp_pos, actual_nonsense]]:
            for position in list(positions[0]):
                print(nonse_chromo, str(int(position) - 1), position, trans, sep='\t', file=positions[1])

    possible_nonsense.close()
    actual_nonsense.close()

if __name__ == '__main__':
    main()
