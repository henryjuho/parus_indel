#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-keep_indels',
                        help='if specified retains indels in alignment',
                        default=False,
                        action='store_true')
    args = parser.parse_args()

    extracted_data = []

    first_line = True
    for line in sys.stdin:

        line = line.split('\t')
        spp = line[4].split(',')
        sequences = line[7].split(',')

        # extracts spp names
        if first_line:
            extracted_data = [[spp_name, ''] for spp_name in spp]
            first_line = False

        # removes indels if desired
        if not args.keep_indels:
            if '-' in ''.join(sequences):
                # need to retain first base as indel tagged on to prev position
                stripped_seqs = ['' for y in sequences]
                for i in range(0, len(sequences[0])):
                    current_seqs = [x[i] for x in sequences]
                    if '-' in current_seqs:
                        continue
                    else:
                        stripped_seqs = [''.join(z) for z in zip(stripped_seqs, current_seqs)]

                sequences = stripped_seqs

        for i in range(0, len(sequences)):
            extracted_data[i][1] += sequences[i]

    # output phylip
    print('\t{}\t{}'.format(len(extracted_data), len(extracted_data[0][1])))
    for i in range(0, len(extracted_data)):
        print(extracted_data[i][0])
        for q in range(0, len(extracted_data[i][1]), 60):
            print(extracted_data[i][1][q: q+60])

if __name__ == '__main__':
    main()
