#!/usr/bin/env python

import argparse
import re
import sys
from scipy import stats

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-data1', help='Results file from full model', required=True)
parser.add_argument('-data2', help='Results file from reduced/null model', required=True)
parser.add_argument('-df', help='Degrees of freedom', required=True, type=int)
args = parser.parse_args()

# variables
mods = {'mod1': args.data1, 'mod0': args.data2}
df = args.df


# functions
def convergence(v1, v2):
    if abs(v1 - v2) / (0.5 * (v1 + v2)) < 10**-5:
        return True
    else:
        return False


def get_variable_columns(header_line):
    column_dict = {header_line[x]: x for x in range(0, len(header_line))}
    return column_dict


def get_var(key, record, cp_dict):
    try:
        return record[cp_dict[key]]
    except KeyError:
        return 'NA'

# process results files
output_lines = []
mod_lnL = {'mod1': 0, 'mod0': 0}
for m in mods.keys():
    with open(mods[m]) as results:
        algorithm = mods[m].split('.')[-2]
        data_lines = [line.rstrip('\n').split('\t') for line in results.readlines() if re.match(r'^\d', line) or
                      line.startswith('run')]
        data = data_lines[1:4]
        header = data_lines[0]

        # column positions
        cp = get_variable_columns(header)

        # check exit code
        if not data[0][2] == data[1][2] == '3':
            sys.exit(m + ' exit code != 3')

        # check convergence
        if convergence(float(data[0][-1]), float(data[1][-1])) is False:
            sys.exit(m + ' failed to converge')

        # extract top result
        snp_theta = get_var('d1_theta', data[0], cp)
        snp_error = get_var('d1_e_1', data[0], cp)
        snp_gamma = get_var('d1_gamma_1', data[0], cp)
        indel_theta = get_var('d2_theta', data[0], cp)
        indel_kappa = get_var('d2_kappa', data[0], cp)
        insertion_gamma = get_var('d2_gamma_1', data[0], cp)
        deletion_gamma = get_var('d2_gamma_2', data[0], cp)
        insertion_error = get_var('d2_e_1', data[0], cp)
        deletion_error = get_var('d2_e_2', data[0], cp)
        indel_theta2 = get_var('d3_theta', data[0], cp)
        indel_kappa2 = get_var('d3_kappa', data[0], cp)
        insertion_gamma2 = get_var('d3_gamma_1', data[0], cp)
        deletion_gamma2 = get_var('d3_gamma_2', data[0], cp)
        insertion_error2 = get_var('d3_e_1', data[0], cp)
        deletion_error2 = get_var('d3_e_2', data[0], cp)

        lnL = get_var('lnL', data[0], cp)

        output_string = '\t'.join([m, algorithm,
                                   snp_theta, snp_error, snp_gamma,
                                   indel_theta, indel_kappa,
                                   insertion_gamma, deletion_gamma,
                                   insertion_error, deletion_error,
                                   indel_theta2, indel_kappa2,
                                   insertion_gamma2, deletion_gamma2,
                                   insertion_error2, deletion_error2,
                                   lnL])

        output_lines.append(output_string)

        # record lnL
        mod_lnL[m] = float(lnL)

# perform ratio test
diff = 2*(mod_lnL['mod1'] - mod_lnL['mod0'])
pval = stats.chi2.sf(diff, df)

# write to standard out
for mod_result in output_lines:
    print(mod_result + '\t' + str(pval))
