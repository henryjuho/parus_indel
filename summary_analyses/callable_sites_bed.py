#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam


def bed_regions(bed_file, chromo):
    open_bed = pysam.TabixFile(bed_file)
    for line in open_bed.fetch(chromo, parser=pysam.asTuple()):
        start = int(line[1])
        end = int(line[2])

        yield start, end


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-call_fa', help='Callable fasta file', required=True)
    parser.add_argument('-chr_list', help='File of chromosomes to calc callable sites for', required=True)
    parser.add_argument('-bed', help='Bed file of regions to count callables sites', required=True)
    parser.add_argument('-tag', help='region name', required=True)
    args = parser.parse_args()

    # variables
    fa = pysam.FastaFile(args.call_fa)
    bed_file = args.bed
    chr_list = [x.rstrip() for x in open(args.chr_list)]
    region = args.tag
    call_data = {'ALL': 0, 'POL': 0}

    # get call sites for all chr and regions
    for chromo in chr_list:
        fasta_string = fa.fetch(chromo)

        if region == 'ALL':
            callable_sites_all = fasta_string.upper().count('K')
            callable_sites_pol = fasta_string.count('K')
        elif region == 'AR':
            callable_sites_all = fasta_string.upper().count('R')
            callable_sites_pol = fasta_string.count('R')

        # handle optional bed files
        else:
            callable_sites_all = 0
            callable_sites_pol = 0
            try:
                for coord_range in bed_regions(bed_file, chromo):
                    fasta_seq = fasta_string[coord_range[0]:coord_range[1]]
                    callable_sites_all += fasta_seq.upper().count('K')
                    callable_sites_pol += fasta_seq.count('K')

            except ValueError:
                pass

        call_data['ALL'] += callable_sites_all
        call_data['POL'] += callable_sites_pol

    # output sites
    print('category', 'variation', 'callable', sep='\t')

    for var in ['SNP', 'INS', 'DEL', 'INDEL']:

        if var == 'SNP' or var == 'INDEL':
            call = call_data['ALL']
        else:
            call = call_data['POL']

        print(args.tag, var, call, sep='\t')


if __name__ == '__main__':
    main()
