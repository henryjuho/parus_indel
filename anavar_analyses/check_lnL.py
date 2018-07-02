from __future__ import print_function
import sys


def main():

    res_files = [x.rstrip() for x in sys.stdin]

    print('rank', 'lnL', 'model')
    for res in res_files:

        counter = 0
        model = res.split('/')[-1].split('.')[0]

        for line in open(res):

            line = line.split()

            if len(line) > 0 and line[0].isdigit():

                counter += 1

                print(counter, line[-1], model)


if __name__ == '__main__':
    main()
