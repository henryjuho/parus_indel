#!/usr/bin/env python

from __future__ import print_function
import argparse
import gzip
import os
import sys
sys.path.insert(0, os.getenv('HOME') + '/parus_indel/annotation')
from degen_to_bed_gt import contig_dict


def main():
    # argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-cds_fa',
                        help='Fasta file of CDS transcript sequences',
                        required=True)
    parser.add_argument('-contig_key',
                        help=argparse.SUPPRESS,
                        default='{}/parus_indel/annotation/contigs_key.txt'.format(os.getenv('HOME')))
    args = parser.parse_args()

    contig_data = contig_dict(args.contig_key)

    skip = False
    prev_chromo = ''
    out_file_stem = args.cds_fa.replace('.fna.gz', '')
    out_file_str = '{}_{}.fna.gz'
    out_file = ''

    for line in gzip.open(args.cds_fa):
        if line.startswith('>'):

            ref_seq_chromo = '_'.join(line.split('[')[0].split('|')[1].split('_')[0:2])

            if ref_seq_chromo not in contig_data.keys():
                skip = True
                continue
            else:
                skip = False

                chromo = contig_data[ref_seq_chromo]

                if prev_chromo != chromo:
                    if type(out_file) != str:
                        out_file.close()

                    chr_out = out_file_str.format(out_file_stem, chromo)
                    out_file = gzip.open(chr_out, 'w')

                prev_chromo = chromo

                print(line, file=out_file)

        else:
            if skip:
                continue
            else:
                print(line, file=out_file)

if __name__ == '__main__':
    main()
