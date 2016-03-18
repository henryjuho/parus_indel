import argparse

# arguments
parser = argparse.ArgumentParser()
parser.add_argument('-fa_list', help='List file of fasta files to merge', required=True)
args = parser.parse_args()

# variables
fastas = [fa.rstrip('\n') for fa in open(args.fa_list)]
output = open('/fastdata/bop15hjb/GT_data/Multispecies_alignment/Whole_genomes_chr_only/'
              'Meleagris_gallopavo.UMD2.dna_sm.fa', 'w')

# merge
for fa in fastas:
    for line in open(fa):
        if line.startswith('>'):
            header = line.split()[0] + '\n'
            print('Processing ' + header)
            output.write(header)
        else:
            output.write(line)

print('Merger complete!')
