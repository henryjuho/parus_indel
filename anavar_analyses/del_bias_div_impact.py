from __future__ import print_function
import anavar_utils as an
import sys


def adjust_div(div_tuple, adjuster, direction):

    if direction == 1:
        adjust = (div_tuple[0]-adjuster, div_tuple[1]+adjuster)
    else:
        adjust = (div_tuple[0] + adjuster, div_tuple[1] - adjuster)

    return adjust


def main():
    res = sys.argv[1]

    results = an.ResultsFile(open(res))

    dn_truth = (0.000121257528807, 0.000200465669246)  # (ins, del)
    ds_truth = (0.00183918200407, 0.00383421025726)  # (ins, del)

    print('percent_bias_adjust', 'alpha', 'var_type', 'rDI_cds', 'rDI_nc', sep='\t')

    # estimate alpha assuming different % errors of div polarisation ie 1% error
    for i in range(0, 26):

        if i == 0:
            dn_adjuster = 0
            ds_adjuster = 0

        else:
            percent_error = i/float(100)

            dn_adjuster = sum(dn_truth) * percent_error
            ds_adjuster = sum(ds_truth) * percent_error

        for direction in [-1, 1]:

            dn_adj = adjust_div(dn_truth, dn_adjuster, direction)
            ds_adj = adjust_div(ds_truth, ds_adjuster, direction)

            ins_alpha = results.get_alpha(dn=dn_adj[0], ds=ds_adj[0], var_type='ins')
            del_alpha = results.get_alpha(dn=dn_adj[1], ds=ds_adj[1], var_type='del')

            print(i * direction, ins_alpha, 'ins', dn_adj[1]/dn_adj[0], ds_adj[1]/ds_adj[0], sep='\t')
            print(i * direction, del_alpha, 'del', dn_adj[1] / dn_adj[0], ds_adj[1] / ds_adj[0], sep='\t')


if __name__ == '__main__':
    main()
