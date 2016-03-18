import argparse
import re

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fa_list', help='List file of fasta files to merge', required=True)
args = parser.parse_args()

# variables
fastas = [fa.rstrip('\n') for fa in open(args.fa_list)]
output = open('/fastdata/bop15hjb/GT_data/Multispecies_alignment/Whole_genomes_chr_only/'
              'FicAlb1.5_genomic.fa', 'w')

# merge
for fa in fastas:
    for line in open(fa):
        if line.startswith('>'):
            try:
                header = '>' + re.search(r'chromosome ([LGEMZ]{0,3}[\d]{0,2}[A]?)', line).group(1) + '\n'
            except:
                header = '>' + re.search(r'linkage group ([LGEMZ]{0,3}[\d]{0,2}[A]?)', line).group(1) + '\n'
            print('Processing ' + header)
            output.write(header)
        else:
            output.write(line)

print('Merger complete!')
