#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse
import subprocess


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-keep_indels',
                        help='if specified retains indels in alignment',
                        default=False,
                        action='store_true')
    parser.add_argument('-out_stem', required=True)
    args = parser.parse_args()

    extracted_data = []
    out_fastas = []

    first_line = True
    for line in sys.stdin:

        line = line.split('\t')
        spp = line[4].split(',')
        sequences = line[7].split(',')

        # extracts spp names
        if first_line:
            fa_names = [args.out_stem + '.' + spp_name + '.fa' for spp_name in spp]
            out_fastas = [open(x, 'w') for x in fa_names]
            for i in range(0, len(spp)):
                print('>{}'.format(spp[i]), file=out_fastas[i])
                extracted_data.append('')
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
            extracted_data[i] += sequences[i]

        # keep on top of printing
        if len(extracted_data[0]) >= 60:
            for i in range(0, len(sequences)):
                print(extracted_data[i][:60], file=out_fastas[i])
                extracted_data[i] = extracted_data[i][60:]

    # print tail seq
    for i in range(0, len(extracted_data)):
        print(extracted_data[i], file=out_fastas[i])

    # close files
    for fa in out_fastas:
        fa.close()

    # cat cmd
    cat = 'cat {z}*.fa > {z}.fa'.format(z=args.out_stem)
    subprocess.call(cat, shell=True)

    # index
    index = 'samtools faidx {}.fa'.format(args.out_stem)
    subprocess.call(index, shell=True)


if __name__ == '__main__':
    main()
