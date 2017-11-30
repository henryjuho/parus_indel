from __future__ import print_function
import sys
import subprocess


vcf = ('/fastdata/bop15hjb/GT_data/BGI_BWA_GATK/Analysis_ready_data/final/'
       'bgi_10birds.filtered_indels.pol.anno.recomb.line.vcf.gz')

bin_files = sorted([(int(x.rstrip().replace('.bed.gz', '').split('.')[-1].replace('bin', '')),
                     x.rstrip()) for x in sys.stdin])

print('bin\tcount\ttype')

for bin_file in bin_files:

    sfs_cmd = ('bedtools intersect -header -a {} -b {} | '
               '~/sfs_utils/vcf2raw_sfs.py -region intergenic -region intron -mode {} | '
               'wc -l')

    for var in ['ins', 'del']:

        n = subprocess.Popen(sfs_cmd.format(vcf, bin_file[1], var),
                             shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[0]

        print(bin_file[0], n, var, sep='\t')

