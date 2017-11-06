#!/usr/bin/env python

from __future__ import print_function
import sys
import re

# make chromosome dict
for line in sys.stdin:
    if line.startswith('#'):
        continue
    elif line.split()[2] == 'region':
        try:
            line = line.split()
            ID = line[0]
            info = '\t'.join(line[8:])
            try:
                chromo = 'chr' + re.search(r';chromosome=([\dABLGEMTWZUnkow]{1,7});', info).group(1)
            except AttributeError:
                try:
                    chromo = 'chr' + re.search(r';Name=([\dABLGEMTWZUnkow]{1,7});', info).group(1)
                except AttributeError:
                    try:
                        chromo = 'chr' + re.search(r';linkage-group=([\dABLGEMTWZUnkowC_]{1,18});', info).group(1)
                    except AttributeError:
                        try:
                            chromo = 'chr' + re.search(r';genome=(mitochondrion);', info).group(1)
                        except AttributeError:
                            continue
            print('\t'.join([chromo, ID]))
        except IndexError:
            continue
