from __future__ import print_function
import sys
import gzip
from vcf2raw_sfs import vcf2sfs

vcf = '/Users/henryjuho/sharc_fastdata/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz'

bin_dict = {}
for bin_file in sys.stdin:
    bin_file = bin_file.rstrip()

    bin_id = bin_file.replace('.bed.gz', '').split('.')[-1].replace('bin', '')
    n_ins = 0
    n_del = 0

    for line in gzip.open(bin_file):
        chromo, start, stop = line.split()[0], int(line.split()[1]), int(line.split()[2])

        ins = vcf2sfs(vcf, mode='ins', chromo=chromo, start=start, stop=stop, regions=['intergenic'], auto_only=True)
        dels = vcf2sfs(vcf, mode='del', chromo=chromo, start=start, stop=stop, regions=['intergenic'], auto_only=True)

        n_ins += len(list(ins))
        n_del += len(list(dels))

    bin_dict[bin_id] = {'ins': n_ins, 'del': n_del}

print('bin\tcount\ttype')
for x in sorted(bin_dict.keys()):
    for z in ['ins', 'del']:
        print(x, bin_dict[x][z], z, sep='\t')
