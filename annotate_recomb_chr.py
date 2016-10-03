#!/usr/bin/env python

import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-vcf', help='/input/directory/vcf', required=True)
parser.add_argument('-ref', help='/path/to/reference/genome.fa', required=True)
parser.add_argument('-poly', help='File containing variable values for polynomial equations for each chromosome',
                    required=True)
parser.add_argument('-chr', help='Chromosome to annotate', required=True)
parser.add_argument('-out', help='/output/directory/', required=True)
args = parser.parse_args()

# variables
vcf = args.vcf
ref_index = args.ref + '.fai'
chr_lengths = {line.split('\t')[0]: float(line.split('\t')[1]) for line in open(ref_index) if line.startswith('chr')}
chromo = args.chr
chromo_length = chr_lengths[chromo]
micro = ['chr' + str(i) for i in range(13, 29) if i != 16 if i != 25]
macro = ['chr' + str(i) for i in range(1, 13)] + ['chrZ']
polynomial_values = {'chr' + x.split()[0]: (float(x.split()[1]), float(x.split()[2]), float(x.split()[3]))
                     for x in open(args.poly)}
try:
    pol_coef = polynomial_values[chromo]
except KeyError:
    pol_coef = 'not_fitted'
out_dir = args.out
annotated_vcf = open(out_dir + vcf[vcf.rfind('/')+1:].replace('.vcf', '.recomb.' + chromo + '.vcf'), 'w')


# functions
def get_mid_25(length):
    distance_from_end = (length * 0.75) / 2.0
    start_point = distance_from_end
    end_point = length - distance_from_end
    return start_point, end_point


def get_end_3mb(length):
    distance_from_end = 3000000
    start_point_1 = 0
    end_point_1 = distance_from_end
    start_point_2 = length - distance_from_end
    end_point_2 = length
    return (start_point_1, end_point_1), (start_point_2, end_point_2)


# Toni's polynomial function to predict recomb
def predict_recomb(pos, x, y, z):
    predicted_value = 3*x*pos**2+2*y*pos+z
    if predicted_value > 0.0:
        return predicted_value
    else:
        return 0.0

# loop through vcf and identify category for each variant
all_variants = 0
previous_line = ''
for line in open(vcf):
    if line.startswith('#'):
        if line.startswith('##contig') and previous_line.startswith('##INFO'):
            new_info = ('##INFO=<ID=RR,Number=1,Type=Float,Description="Annotation of recombination rate">\n'
                        '##INFO=<ID=RBIN,Number=1,Type=String,Description="Annotation of crude recombination bin">\n')
            annotated_vcf.write(new_info)
        previous_line = line
        annotated_vcf.write(line)
    elif line.split('\t')[0] == chromo:
        all_variants += 1
        split_line = line.split('\t')
        chromo, start, info = split_line[0], int(split_line[1])-1, split_line[7]
        if pol_coef == 'not_fitted':
            recomb_rate = 'None'
        else:
            recomb_rate = predict_recomb(start, pol_coef[0], pol_coef[1], pol_coef[2])
        if recomb_rate != 'None':
            info += ';RR=' + str(recomb_rate)
        recomb_bin = 'None'
        if chromo in micro:
            recomb_bin = 'micro'
        elif chromo in macro:
            mid_25 = get_mid_25(chromo_length)
            ends = get_end_3mb(chromo_length)
            if mid_25[0] <= start <= mid_25[1]:
                recomb_bin = 'macro_mid'
            elif ends[0][0] <= start <= ends[0][1] or ends[1][0] <= start <= ends[1][1]:
                recomb_bin = 'macro_end'
        if recomb_bin != 'None':
            info += ';RBIN=' + recomb_bin

        # write annotated line
        annotated_vcf.write('\t'.join(split_line[0:7]) + '\t' + info + '\t' + '\t'.join(split_line[8:]))

    else:
        continue

annotated_vcf.close()
