#!/usr/bin/env python

from __future__ import print_function
import argparse
from pysam import VariantFile
import pysam
import gzip


def polarisable(vcf_var, wga_bed):
    chrom, pos, ref, alt = vcf_var.contig, vcf_var.pos, vcf_var.ref, vcf_var.alts

    if alt is None:
        var_type = 'MONO'
        # catch monomorphic sites recorded as INDELs in vcf
        ref = ref[0]

    elif len(ref) > 1:
        var_type = 'INDEL'

    else:
        var_type = 'SNP'

    # get aligned seq for var
    try:
        var_align = [x for x in wga_bed.fetch(chrom, pos - 1, pos - 1 + len(ref),
                                              parser=pysam.asTuple())]

    # catch if whole chromo not in align
    except ValueError:
        return False, 'chr_missing'

    # skip if not in alignment
    if len(var_align) == 0:
        return False, 'not_aligned'

    elif len(var_align) > 1:
        if var_type != 'INDEL':
            return False, 'in_indel'

        # catch deletions rel to ref that uniq to ref spp
        else:
            # merge sequences from multiple bed rows
            merged_align = [''.join(y) for y in zip(*[x[7].split(',') for x in var_align])]
            if '-' not in ''.join(merged_align):
                var_align = merged_align
            else:
                return False, 'indel_hotspot'

    else:
        var_align = [x[7] for x in var_align][0].split(',')

    var_align = [x.upper() for x in var_align]
    # print chrom, pos, len(ref), ref, alt, var_align

    # skip positions without full coverage
    if '?' in ''.join(var_align):
        return False, 'low_coverage'

    # skip indel hotspots
    indel_sequences = [y.rstrip('-') for y in var_align]
    if '-' in ''.join(indel_sequences):
        return False, 'indel_hotspot'

    # skips sites where ref allele differs from that in alignment, ie insertion within INDEL
    if ref != var_align[0].upper():
        return False, 'indel_hotspot'

    else:
        # identify if ref or alt is ancestral
        out_group_seqs = var_align[1:]

        # deal with monomorphic sites
        if var_type == 'MONO':
            if len(set(out_group_seqs)) == 1:
                return True, 'polarisable_mono'
            else:
                return False, 'ambiguous_mono'

        ref_anc = True
        alt_anc = True

        # identify ref ancestral
        for sequence in out_group_seqs:
            if var_type == 'INDEL':
                if len(ref) != len(sequence.rstrip('-')):
                    ref_anc = False
                    break
            else:
                if ref != sequence.rstrip('-'):
                    ref_anc = False
                    break

        # identify alt ancestral
        for sequence in out_group_seqs:
            if var_type == 'INDEL':
                if len(alt[0]) != len(sequence.rstrip('-')):
                    alt_anc = False
                    break
            else:
                if alt[0] != sequence.rstrip('-'):
                    alt_anc = False
                    break

        # skip ambiguous sites
        if alt_anc is ref_anc:
            return False, 'ambiguous'

        # can be polarised
        return True, 'polarisable'


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf',
                        help='Allsites vcf to apply filters to and get callable sites',
                        required=True)
    parser.add_argument('-bed', '--bed_repeats',
                        help='BED file with repeat regions listed',
                        required=True)
    parser.add_argument('-ar_bed', '--ar_bed',
                        help='BED file of ancestral repeats',
                        default='None')
    parser.add_argument('-DF', '--DepthFilter',
                        help='Defines abnormal depth eg) 2 means abnormal depth is twice and half the mean depth',
                        default=2.0,
                        type=float)
    parser.add_argument('-mean_depth', '--mean_depth',
                        help='Mean coverage depth of samples',
                        default=44.0)
    parser.add_argument('-N', '--no_individuals',
                        help='Number of individuals in allsites VCF',
                        type=float,
                        default=10.0)
    parser.add_argument('-chr',
                        help='Specifies chromosome to extract callable sites for, if ALL will run a job for each, '
                             '-chr ALL can only be specified in conjunction with -sub',
                        default='ALL')
    parser.add_argument('-start',
                        help='start coordinates to get callable sites for',
                        type=int,
                        required=True)
    parser.add_argument('-end',
                        help='end coordinate for region',
                        type=int,
                        required=True)
    parser.add_argument('-pol',
                        help='If specified will check if site can be polarised, takes a wga bed file',
                        default='None')
    args = parser.parse_args()

    # variables
    all_sites = args.vcf
    repeat_bed = args.bed_repeats
    line_bed = args.ar_bed
    start = args.start
    stop = args.end
    filter_factor = args.DepthFilter
    all_data_mean_depth = float(args.mean_depth)
    no_indiv = args.no_individuals
    chromosome = args.chr
    pol = args.pol

    # calculate depth cutoffs
    lower_depth_limit = all_data_mean_depth / filter_factor
    upper_depth_limit = all_data_mean_depth * filter_factor

    repeats = set()
    # get bed regions per chromo
    for x in open(repeat_bed):
        if x.split()[0] == chromosome:
            repeats |= {y for y in range(int(x.split()[1]), int(x.split()[2]))}

    lines = set()
    # get bed regions per chromo
    if line_bed != 'None':
        for x in gzip.open(line_bed):
            if x.split()[0] == chromosome:
                lines |= {y for y in range(int(x.split()[1]), int(x.split()[2]))}

    # loop through allsites for chromosome
    counter = 0
    fasta_string = '>{}:{}-{}'.format(chromosome, start, stop)
    if pol != 'None':
        wga_bed = pysam.TabixFile(pol)
    else:
        wga_bed = None

    # move through vcf
    print(fasta_string)
    fasta_string = ''
    prev_position = 0
    for line in VariantFile(all_sites).fetch(chromosome, start+1, stop):

        # catch missing sites in allsites (new gatk3.7 feature)
        position = int(line.pos)
        diff = position - prev_position
        if diff != 1:
            missed_bases = ''.join(['1' for i in range(0, diff-1)])
            fasta_string += missed_bases
        prev_position = position

        # add line break every 60 bases
        if len(fasta_string) >= 60:
            if len(fasta_string) == 60:
                print(fasta_string)
                fasta_string = ''
            else:
                print(fasta_string[:60])
                fasta_string = fasta_string[60:]
        counter += 1

        # check for ns
        if line.ref == 'N':
            fasta_string += '0'
            continue

        # depth filter
        try:
            cumulative_depth = line.info["DP"]
        except KeyError:
            fasta_string += '1'
            continue

        locus_mean_depth = cumulative_depth / no_indiv
        if lower_depth_limit <= locus_mean_depth <= upper_depth_limit:

            # repeat filter
            if line.pos not in repeats:

                # check if polarisable
                if pol != 'None':
                    can_polarise = polarisable(line, wga_bed)[0]
                    if can_polarise is False:
                        fasta_string += 'k'
                        continue
                    else:
                        fasta_string += 'K'
                        continue
                else:
                    fasta_string += 'k'
                    continue

            else:
                if line.pos in lines:

                    # check if polarisable
                    if pol != 'None':
                        can_polarise = polarisable(line, wga_bed)[0]
                        if can_polarise is False:
                            fasta_string += 'r'
                            continue
                        else:
                            fasta_string += 'R'
                            continue
                    else:
                        fasta_string += 'r'
                        continue
                else:
                    fasta_string += '1'
                continue

        else:
            fasta_string += '1'
            continue

    print(fasta_string)


if __name__ == '__main__':
    main()
