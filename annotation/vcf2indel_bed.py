from __future__ import print_function
import sys
import argparse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-read_frame', help='Consider reading frame', choices=['inframe', 'shift'], default='ALL')
    args = parser.parse_args()

    for line in sys.stdin:

        chromo, pos, ref_len, alt_len = line.split()[0], int(line.split()[1]), \
                                        len(line.split()[3]), len(line.split()[4])

        indel_len = abs(ref_len - alt_len)

        if args.read_frame == 'inframe' and indel_len % 3 != 0:
            continue

        elif args.read_frame == 'shift' and indel_len % 3 == 0:
            continue

        print(chromo, pos-1, pos-1+ref_len, sep='\t')


if __name__ == '__main__':
    main()
