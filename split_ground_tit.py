#!/usr/bin/env python

ground_tit_fasta = '/fastdata/bop15hjb/GT_data/Multispecies_alignment/Alignment_genomes/Groundtit/' \
                   'GCF_000331425.1_PseHum1.0_genomic.fna'
out_dir = '/fastdata/bop15hjb/GT_data/Multispecies_alignment/Alignment_genomes/Groundtit/'

# count scaffolds
no_scaf = 0.0

with open(ground_tit_fasta) as fasta_in:
    for line in fasta_in:
        if line.startswith('>'):
            no_scaf += 1.0

scaffolds_per_bin = (no_scaf / 300.0) + 1

# write binned scaffold fastas
output_tracker = 1
scaffold_tracker = 0.0
output_prefix = 'Groundtit_scaffold_bin_'
out_locations = [[out_dir + output_prefix + str(output_tracker) + '.fa', output_prefix + str(output_tracker)]]
output_fasta = open(out_dir + output_prefix + str(output_tracker) + '.fa', 'w')
with open(ground_tit_fasta) as fasta_in:
    for line in fasta_in:
        if line.startswith('>'):
            if scaffold_tracker > scaffolds_per_bin:
                scaffold_tracker = 0.0
                output_tracker += 1
                output_fasta.close()
                output_fasta = open(out_dir + output_prefix + str(output_tracker) + '.fa', 'w')
                out_locations.append([out_dir + output_prefix + str(output_tracker) + '.fa',
                                      output_prefix + str(output_tracker)])
            scaffold_tracker += 1.0
        output_fasta.write(line)

output_fasta.close()

# write chromo_list
with open(out_dir + 'Groundtit_scaffold_bins.list', 'w') as chromo_list:
    for bins in out_locations:
        chromo_list.write(bins[0] + '\t' + bins[1] + '\n')

print('Done')
