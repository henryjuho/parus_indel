#!/usr/bin/env python

from __future__ import print_function
import sys


def fasta_print(contig, seq):

    print('>{}'.format(contig))

    for i in range(0, len(seq), 60):
        print(seq[i:i+60])


def main():

    sequence = ''
    chromo = ''
    prev_chromo = ''

    for file_name in sys.stdin:
        for line in open(file_name.rstrip()):
            if line.startswith('>'):
                prev_chromo = chromo
                chromo = line.split(':')[0].replace('>', '')
                if prev_chromo == chromo or prev_chromo == '':
                    continue
                else:
                    fasta_print(prev_chromo, sequence)
                    sequence = ''
            else:
                sequence += line.rstrip()

    fasta_print(prev_chromo, sequence)

if __name__ == '__main__':
    main()
