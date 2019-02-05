from __future__ import print_function
import sys


def main():

    res_files = [x.rstrip() for x in sys.stdin]

    print('rank', 'lnL', 'model', sep='\t')
    for res in res_files:

        lnl_list = []
        model = res.split('/')[-1].split('.')[0]

        for line in open(res):

            line = line.split()

            if len(line) > 0 and line[0].isdigit():

                if float(line[-1]) > 0:

                    lnl_list.append(line[-1])

                else:

                    lnl_list.append('NA')

        for i in range(1, len(lnl_list)+1):

            print(i, lnl_list[-i], model, sep='\t')


if __name__ == '__main__':
    main()
