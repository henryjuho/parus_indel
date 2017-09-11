#!/usr/bin/env python

from __future__ import print_function
import argparse
import pysam
import itertools


def ranges(i):
    for a, b in itertools.groupby(enumerate(i), lambda (x, y): y - x):
        b = list(b)
        yield b[0][1], b[-1][1]


def bed_regions(bed_file, chromo):
    open_bed = pysam.TabixFile(bed_file)
    coords = []
    for line in open_bed.fetch(chromo, parser=pysam.asTuple()):
        start = int(line[1])
        end = int(line[2])
        coords += range(start, end)

    return coords


def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-call_fa', help='Callable fasta file', required=True)
    parser.add_argument('-chr_list', help='File of chromosomes to calc callable sites for', required=True)
    parser.add_argument('-opt_bed', help='Optional bed file of regions to count callables sites, with associated label'
                                         'i.e. /path/to/file.bed.gz,my_sub_region', action='append')
    args = parser.parse_args()

    # variables
    fa = pysam.FastaFile(args.call_fa)
    bed_files = {x.split(',')[1]: x.split(',')[0] for x in args.opt_bed}
    chr_list = [x.rstrip() for x in open(args.chr_list)]
    regions = ['ALL', 'AR'] + bed_files.keys()
    call_data = {x: {y: {'ALL': 0, 'POL': 0} for y in regions} for x in ['ALL'] + chr_list}

    # {chromo: {all: 0, CDS: 0, intron: 0 ...}}

    # get call sites for all chr and regions
    for chromo in chr_list:
        fasta_string = fa.fetch(chromo)
        for region in regions:
            if region == 'ALL':
                callable_sites_all = fasta_string.upper().count('K')
                callable_sites_pol = fasta_string.count('K')
            elif region == 'AR':
                callable_sites_all = fasta_string.upper().count('R')
                callable_sites_pol = fasta_string.count('R')

            # handle optional bed files
            else:
                degen_bed = bed_files[region]
                try:
                    degen_coords = bed_regions(degen_bed, chromo)

                    callable_sites_all = ''.join([fasta_string[x] for x in degen_coords]).upper().count('K')
                    callable_sites_pol = ''.join([fasta_string[x] for x in degen_coords]).count('K')
                except ValueError:
                    callable_sites_all = 0
                    callable_sites_pol = 0

            call_data[chromo][region]['ALL'] += callable_sites_all
            call_data[chromo][region]['POL'] += callable_sites_pol

            call_data['ALL'][region]['ALL'] += callable_sites_all
            call_data['ALL'][region]['POL'] += callable_sites_pol

    # output sites
    print(','.join(['contig', 'region', 'all_callable', 'pol_callable']))
    for seq in sorted(call_data.keys()):
        for reg in call_data[seq].keys():
            print(','.join([seq, reg, str(call_data[seq][reg]['ALL']), str(call_data[seq][reg]['POL'])]))

if __name__ == '__main__':
    main()
