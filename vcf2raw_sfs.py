#!/usr/bin/env python

import argparse
import pysam
import sys


# functions
def get_derived_freq(vcf_line, run_mode):

    """
    function takes a pysam vcf variant and a variant mode and returns the derived allele frequency for that variant,
    if it matches the given variant mode (eg ins, if insertion spectrum required)
    :param vcf_line: pysam variant
    :param run_mode: string
    :return: None or float
    """

    ref_seq = vcf_line.ref
    alt_seq = vcf_line.alts[0]
    alt_freq = round(vcf_line.info['AF'][0], 2)
    try:
        anc_seq = vcf_line.info['AA']
    except KeyError:  # ie. not polarised
        return None

    # set derived sequence and freq
    if alt_seq == anc_seq:
        derv_seq = ref_seq
        derv_freq = 1 - alt_freq
    elif ref_seq == anc_seq:
        derv_seq = alt_seq
        derv_freq = alt_freq
    else:
        return None

    # determine type
    if len(anc_seq) > len(derv_seq):  # deletion
        if run_mode == 'del':
            return derv_freq
        else:
            return None
    elif len(anc_seq) < len(derv_seq):  # insertion
        if run_mode == 'ins':
            return derv_freq
        else:
            return None
    else:  # snp
        return derv_freq


def get_minor_freq(vcf_line):

    """
    takes a pysam variant and returns the minor allele frequency
    :param vcf_line: pysam variant
    :return: float
    """

    alt_allele_freq = round(vcf_line.info['AF'][0], 2)
    if alt_allele_freq <= 0.5:
        return alt_allele_freq
    else:
        return 1 - alt_allele_freq


def get_out_freq(vcf_line, pol, run_mode):

    """
    takes pysam variant, polarisation argument and variant run type
    :param vcf_line: pysam variant
    :param pol: bool
    :param run_mode: str
    :return: None or float
    """

    if pol is True:
        return get_derived_freq(vcf_line, run_mode)
    else:
        return get_minor_freq(vcf_line)


def in_regions(vcf_line, target_regions):

    """
    takes a pysam variant and sees if it falls within a specified genomic region
    :param vcf_line: pysam variant
    :param target_regions: list or None
    :return: bool
    """

    if target_regions is None:
        return True
    else:
        try:
            var_region = vcf_line.info['ANNO']
            if var_region in target_regions:
                return True
            else:
                return False
        except KeyError:
            return False


def is_degen(vcf_line, target_degen):

    """
    takes a pysam variant and desired degeneracy and returns true or false for that variant
    :param vcf_line: pysam variant
    :param target_degen: int
    :return: bool
    """

    if target_degen is None:
        return True
    else:
        try:
            degeneracy = int(vcf_line.info['DEGEN'])
            if target_degen == degeneracy:
                return True
            else:
                return False
        except KeyError:
            return False


def is_mute_type(vcf_line, mute_list, pol):

    """
    takes pysam variant and determins if variant is of type listed in mutation list
    :param vcf_line: pysam variant
    :param mute_list: list
    :param pol: bool
    :return:
    """
    # strong = CG weak = AT
    base_types = {'A': 'W', 'T': 'W', 'C': 'S', 'G': 'S'}
    if mute_list is None:
        return True
    else:
        ref_base = base_types[vcf_line.ref]
        alt_base = base_types[vcf_line.alts[0]]
        if ref_base == alt_base == 'W':
            mutation_type = 'WW'
        elif ref_base == alt_base == 'S':
            mutation_type = 'SS'
        else:
            if pol is False:
                return False
            else:
                try:
                    anc_base = base_types[vcf_line.info['AA']]
                    if anc_base == ref_base:
                        mutation_type = anc_base + alt_base
                    else:
                        mutation_type = alt_base + ref_base
                except KeyError:
                    return False
        if mutation_type in mute_list:
            return True
        else:
            return False


def allele_num_ok(vcf_line, no_samples, multi):

    """
    checks to see if variant is biallelic
    :param vcf_line: pysam variant
    :param no_samples: int
    :param multi: bool
    :return: bool
    """

    if multi is False:
        pos_biallelic_freqs = [i/float(2*no_samples) for i in range(1, 2*no_samples)]
        alt_allele_freq = round(vcf_line.info['AF'][0], 2)
        if alt_allele_freq in pos_biallelic_freqs:
            return True
        else:
            return False
    else:
        return True


# main call
def main():

    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-vcf', help='VCF file to extract sfs from', required=True)
    parser.add_argument('-chr', help='Chromosome to extract', default='ALL')
    parser.add_argument('-region', help='Genomic regions to extract, default = ALL', action='append')
    parser.add_argument('-mode', help='Variant mode to run in', choices=['snp', 'ins', 'del', 'indel'], required=True)
    parser.add_argument('-degen', help='Degeneracy of coding SNPs to extract (must run with -mode snp',
                        choices=[0, 2, 3, 4], type=int)
    parser.add_argument('-mute_type', help='Mutation type, use only with mode -snp',
                        choices=['WW', 'SS', 'SW', 'WS'], action='append')
    parser.add_argument('-folded', help='If specified will output minor allele spectrum',
                        default=False, action='store_true')
    parser.add_argument('-multi_allelic', help='If specified will not restrict output to biallelic sites',
                        default=False, action='store_true')
    args = parser.parse_args()

    # variables
    vcf_file = pysam.VariantFile(args.vcf)
    chromo = args.chr
    regions = args.region
    fold = args.folded
    mode = args.mode
    degen = args.degen
    mute_type = args.mute_type
    multi_allelic = args.multi_allelic

    # check commandline options
    if mode == 'indel' and fold is False:
        sys.exit('-mode indel must be run in conjunction with -folded')
    if mode == 'ins' and fold is True or mode == 'del' and fold is True:
        sys.exit('-mode ins and -mode del cannot be run in conjunction with -folded')
    if degen is not None and mode != 'snp':
        sys.exit('-degen can only be specified in conjunction with -mode snp')
    if mute_type is not None and mode != 'snp':
        sys.exit('-mute_type can only be run with -mode snp')

    # loop through vcf
    if chromo == 'ALL':
        vcf = vcf_file.fetch()
    else:
        vcf = vcf_file.fetch(chromo)

    number_samples = len(vcf_file.header.samples)

    for variant in vcf:

        # gets relevant freq, minor or derived, see functions
        frequency = get_out_freq(variant, not fold, mode)
        if frequency is None:  # skips when no freq returned, ie unpolarised or wrong var type
            continue

        # gets variant region if regional sfs required
        falls_in_regions = in_regions(variant, regions)

        # gets degeneracy if required
        degen_ok = is_degen(variant, degen)

        # gets mutation type if required
        mutetype_ok = is_mute_type(variant, mute_type, not fold)

        # checks if is biallelic
        alleles_ok = allele_num_ok(variant, number_samples, multi_allelic)

        # outputs if all criteria ok
        if falls_in_regions is True and degen_ok is True and mutetype_ok is True and alleles_ok is True:
            print frequency


if __name__ == '__main__':
    main()
