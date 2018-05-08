#!/usr/bin/env python

import argparse
from qsub import *


def file_location(file_string):
    """
    function that strip the filename off a path

    :param file_string: str
    :return: path str
    """
    return file_string[:file_string.rfind('/')+1]


def extract_filename(file_string):
    """
    function that strips a path of a filename

    :param file_string: str
    :return: file str
    """
    return file_string[file_string.rfind('/')+1:]


def main():
    # arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-maf', help='MAF file containg alignment, must be compressed', required=True)
    parser.add_argument('-ref_sp', help='Species to use for coordinates in output bed', required=True)
    parser.add_argument('-ref_fa', help='Fasta file to extract chromosomes for ref species', required=True)
    parser.add_argument('-out_pre', help='Output path and filename prefix', required=True)
    parser.add_argument('-evolgen', help='If specified will run on lab queue', default=False, action='store_true')
    args = parser.parse_args()

    # variables
    maf = args.maf
    ref_sp = args.ref_sp
    chr_grep = 'grep "^>" ' + args.ref_fa
    chromo_list = [x.strip('>') for x in
                   subprocess.Popen(chr_grep, shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')
                   if x.startswith('>chr') and 'random' not in x]
    maf_stem = args.out_pre
    evolgen = args.evolgen

    # submit per chromosome jobs
    out_list_ordered = []
    for chromo in chromo_list:
        chromo_bed_out = maf_stem + '.' + chromo + '.wga.bed.gz'
        out_list_ordered.append(chromo_bed_out)
        maf_to_bed_cmd = ('~/WGAbed/maf_to_bed.py '
                          '-i ' + maf + ' '
                          '-r ' + ref_sp + ' '
                          '-c ' + chromo + ' | '
                          'sort -T ' + file_location(maf) + ' -k1,1 -k2,2n | '
                          'bgzip -c > ' + chromo_bed_out)
        q_sub([maf_to_bed_cmd], out=maf_stem + '.' + chromo,
              jid=ref_sp + '_' + chromo + '.sh', t=48, evolgen=evolgen)

    # cat into whole genome
    cat_cmd = 'zcat ' + ' '.join(out_list_ordered) + ' | bgzip -c > ' + maf_stem + '.wga.bed.gz'
    q_sub([cat_cmd], out=maf_stem + '.cat', hold=[ref_sp + '_' + z + '.sh' for z in chromo_list], evolgen=evolgen)

if __name__ == '__main__':
    main()
