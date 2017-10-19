from __future__ import print_function
import sys


"""
converts variant annotations in GT LINE vcf to format to work with drosophila scripts
"""

for line in sys.stdin:
    line = line.rstrip()
    if line.startswith('#'):
        print(line)
    elif 'ANNO=intergenic' in line:
        print(line.replace('ANNO=intergenic', 'ANNO=intergenic_ar'))
    elif 'ANNO=intron' in line:
        print(line.replace('ANNO=intron', 'ANNO=intron_ar'))
    else:
        print(line)
